import json
import functools
from dataclasses import dataclass
from abc import ABC, abstractmethod
# import cProfile

# Requirements
import numpy as np
from pyxyz import Confpool  # use `pip install pyxyz24 on Linux`
import networkx as nx
from networkx.algorithms import isomorphism
from skspatial.objects import Line, Plane
import scipy
import tqdm

from typing import Union

try:
    from numba import jit
except ImportError:
    print(
        "Calculation will be slower without Numba. Install it with 'pip install numba'"
    )

    def jit(func):

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            return func(*args, **kwargs)

        return wrapper

# STARTING_XYZ = r'C:\Users\artem\OneDrive\_Work_Kunitsyn\Projects\PolyBMSTU\visuals\visual3d\PES-185\PES-185-Wn\PES-185-8d_REMD250-Wn.xyz'
# STARTING_XYZ = r'C:\Users\artem\OneDrive\_Work_Kunitsyn\Projects\PolyBMSTU\visuals\visual3d\PES-185\PES-185-Wn\uniform_polymer.xyz'
# STARTING_DATAFILE = r'C:\Users\artem\OneDrive\_Work_Kunitsyn\Projects\PolyBMSTU\visuals\visual3d\PES-185\PES-185-Wn\uniform_polymer_soft.lmps'

# STARTING_XYZ = 'test/structures.xyz'
# STARTING_DATAFILE = 'test/uniform_polymer_soft.lmps'

STARTING_XYZ = '/home/md/md/PES-210_dens/16.07-250_2500-Mn-T600/alpha_0.5_iter1_restart/traj.xyz'
STARTING_DATAFILE = '/home/md/md/PES-210_dens/16.07-250_2500-Mn-T600/uniform_polymer_soft.lmps'
INDEX = 2500



BOND_LENGTH_CUTOFF = 5.0  # Angstroms


@dataclass
class CellParameters:
    x_vector: np.ndarray
    y_vector: np.ndarray
    z_vector: np.ndarray
    base_vertex: np.ndarray

    _magic_x: Union[float,None] = None
    _magic_y: Union[float,None] = None
    _magic_z: Union[float,None] = None

    def __post_init__(self):
        self._magic_x = float(np.dot(self.x_vector, self.x_vector) / 2)
        self._magic_y = float(np.dot(self.y_vector, self.y_vector) / 2)
        self._magic_z = float(np.dot(self.z_vector, self.z_vector) / 2)


class LAMMPSLoader(ABC):

    def __init__(
        self,
        datafile: str,
        show_progress: bool = True,
        debug_plotting: bool = False,
    ) -> None:
        """Initialize from datafile. Bonds can be parsed at this point

        Args:
            datafile (str): File from LAMMPS (I guess). It should contain sections 'Bonds' (see examples).
            debug_plotting (bool, optional): If debug plotting is enabled, an additional JSON file can be generated. This file will contain all key geometric step that can be visualized in Blender (see script `blender_plot.py`). Defaults to False.
        """
        self.bonds = self.read_bonds(datafile)
        self.datafile = datafile
        self.start_cell_data: CellParameters = self.read_starting_cell(
            datafile)
        self.symbols = None
        self.molgraph = None
        self.show_progress = show_progress

        self.molgraph = self.create_molgraph()
        self.graph_preprocessing()

        self.debug_plotting = debug_plotting
        if debug_plotting:
            self.plot_data = {}

    def read_bonds(self, filename: str) -> list[tuple[int, int]]:
        """Reads topology (i.e. 'Bonds' section) of LAMMPS input file (Right?)

        Args:
            filename (str): Path to the file

        Returns:
            list[tuple[int, int]]: list of pairs of atoms (indexing starts from 0)
        """
        with open(filename, 'r') as f:
            lines = f.readlines()

        bonds_section_title = 'Bonds\n'
        assert bonds_section_title in lines
        bonds_section_index = lines.index(bonds_section_title)

        bonds = []
        get_atom_label = lambda label: int(label) - 1
        for i in range(bonds_section_index + 2, len(lines)):
            parts = lines[i].split()
            if len(parts) != 4:
                break
            bonds.append((
                get_atom_label(parts[2]),
                get_atom_label(parts[3]),
            ))
        return bonds

    def read_starting_cell(
        self,
        filename: str,
    ) -> CellParameters:
        """_summary_
        LAMMPS file header example (with line numbers):
        ```
        1  pySIMM System Object
        15 -11.882100 11.882100 xlo xhi
        16 -11.882100 11.882100 ylo yhi
        17 -11.882100 11.882100 zlo zhi
        ```

        Args:
            datafile (str): _description_

        Returns:
            CellParameters: _description_
        """

        with open(filename, 'r') as f:
            lines = f.readlines()

        def extract_vector(
            dim: int,
            line: str,
            dim_str: str,
        ) -> tuple[np.ndarray, float]:
            # Verify the Nlo Nhi ending of the string
            assert f"{dim_str}lo {dim_str}hi" in line

            # Meaning of hi and lo: https://docs.lammps.org/Howto_triclinic.html
            parts = line.split()
            lo = float(parts[0])
            hi = float(parts[1])

            direction = np.zeros(3)
            # dim = 0, 1 or 2
            direction[dim] = hi - lo
            # (edge of the cell, coordinate of the cell vertex)
            return direction, lo

        x_vector, x_vertex = extract_vector(0, lines[14], 'x')
        y_vector, y_vertex = extract_vector(1, lines[15], 'y')
        z_vector, z_vertex = extract_vector(2, lines[16], 'z')
        base_vertex = np.array([x_vertex, y_vertex, z_vertex])
        return CellParameters(x_vector, y_vector, z_vector, base_vertex)

    def create_molgraph(self) -> nx.Graph:
        """Construct molecular graph using bonds taken topology of LAMMPS input file.
        WARNING: the graph has no element symbols at this point!
        """
        graph = nx.Graph()
        graph.add_edges_from(self.bonds)
        return graph

    @abstractmethod
    def graph_preprocessing(self):
        """`CrossingDetector` and `MoleculeReconstructor` require very different
        graph preprocessing. So this method is abstract
        """
        ...

    def crosscheck_symbols(self, symbols: list[str]) -> None:
        assert self.molgraph is not None
        if self.symbols is None:
            self.symbols = symbols
            for i, sym in enumerate(symbols):
                assert self.molgraph.has_node(i)
                assert 'sym' not in self.molgraph.nodes[i]
                self.molgraph.nodes[i]['sym'] = sym
        else:
            assert len(symbols) == self.molgraph.number_of_nodes()
            current_symbols = [
                self.molgraph.nodes[i]['sym']
                for i in range(self.molgraph.number_of_nodes())
            ]
            assert current_symbols == symbols, f"current_symbols={current_symbols}, given_symbols={symbols}"

    @staticmethod
    @jit
    def mod_vector(p: np.ndarray, q: np.ndarray, mq: float):
        while np.dot(p, q) > mq:
            p -= q
        while np.dot(p, q) < -mq:
            p += q

    @staticmethod
    def mod_cell(p: np.ndarray, cell: CellParameters) -> np.ndarray:
        """Modifies the 'p' object itself!
        The input 'p' and the returned vector are the same object
        """
        LAMMPSLoader.mod_vector(p, cell.x_vector, cell._magic_x)
        LAMMPSLoader.mod_vector(p, cell.y_vector, cell._magic_y)
        LAMMPSLoader.mod_vector(p, cell.z_vector, cell._magic_z)
        return p

    @staticmethod
    def periodic_distance(a: np.ndarray, b: np.ndarray, cell: CellParameters):
        diff_vector = a - b
        LAMMPSLoader.mod_cell(diff_vector, cell)
        return np.linalg.norm(diff_vector)


class CrossingDetector(LAMMPSLoader):

    def graph_preprocessing(self) -> None:
        print('Doing topology analysis...')
        self.ring_indices = self.find_rings()
        print(f'Done topology analysis. Found {len(self.ring_indices)} rings')

    def find_rings(self) -> list[tuple[int, ...]]:
        """Identify rings in molecular graph using a general ring-finding algorithm.

        Returns:
            list[tuple[int, ...]]: each element in a list of atom indices of a unique ring. Indexing starts from 0.
        """

        assert self.molgraph is not None
        bridges = set(nx.bridges(self.molgraph))
        nonbridges = [
            edge for edge in self.molgraph.edges if edge not in bridges
        ]
        cyclic_nodes = {node for edge in nonbridges for node in edge}
        cyclic_graph: nx.Graph = self.molgraph.subgraph(cyclic_nodes)
        assert cyclic_graph.number_of_nodes() > 0, (
            "No rings found in the polymer")

        components = [
            cyclic_graph.subgraph(c)
            for c in nx.connected_components(cyclic_graph)
        ]
        rings = [
            ring for component in components
            for ring in nx.minimum_cycle_basis(component)
        ]

        # This check passes (i.e. ring atoms are in the correct order)
        for ring in rings:
            for bond in zip(ring, ring[1:]):
                assert self.molgraph.has_edge(*bond)
            assert self.molgraph.has_edge(ring[0], ring[-1])
        return rings

    def analyze_xyz_file(
        self,
        xyz_filename: str,
        index = None,
        alpha = None,
        xyz = None,
        sym = None,
    ) -> list[list[dict[str, tuple[int, ...]]]]:
        """This core method of the class. Finds intersections between bonds and rings for each structure in XYZ file.
        Structures and atom indexing in XYZ file must agree with indexing in datafile provided in initialization.
        The structure of the output is as follows:
        {Index of conformation in XYZ file} => List of intersections, where each intersection is represented with a dict
            {
                'bond': tuple(atomA, atomB), # Indexing starts from 1
                'ring': tuple of ring atoms, # Indexing starts from 1
            }

        Args:
            xyz_filename (str): Path to XYZ file to be analyzed.

        Returns:
            list[list[dict[str, tuple[int, ...]]]]: List of intersections descriptions (see above)
        """
        
        _index=index
        p = Confpool()
        if xyz is None:
            p.include_from_file(xyz_filename)
        else:
            p.include_from_xyz(xyz,'--')
            p.atom_symbols = sym
            index=0
            

        result = []
        num_conformers = len(p)

        if index==None:
            if self.show_progress:
                pool_iter = enumerate(tqdm.tqdm(p))
            else:
                pool_iter = enumerate(p)
            for i, m in pool_iter:
                if self.debug_plotting:
                    self.scene_id = i
                    self.temp_intersections = []

                if not self.show_progress:
                    print(f"Processing conformation {i} out of {num_conformers}")
                try:
                    result.append(
                        self.find_conflicts(
                            coords=m.xyz,
                            symbols=p.atom_symbols,
                            # For now, we test on fixed cell parameters
                            cell=self.start_cell_data,
                        ))
                except:
                    print(f'exception at frame {_index}, alpha {alpha}')
                    return None,_index,alpha
                # break

                if self.debug_plotting:
                    from visualization_pipelines.blender_scene import (
                        SimplePoint,
                        Cylinder,
                        ArrowMaterial,
                    )

                    cell = self.start_cell_data
                    cell_vertex_mat = ArrowMaterial(color='#FF0000',
                                                    name='cell_vertex')
                    cell_edge_mat = ArrowMaterial(color='#FFFF00',
                                                name='cell_edge')
                    self.plot_data[i] = [
                        *[
                            SimplePoint(cell.base_vertex + shift,
                                        material=cell_vertex_mat,
                                        size=0.5)
                            for shift in (
                                np.zeros(3), cell.x_vector, cell.y_vector,
                                cell.z_vector, cell.x_vector + cell.y_vector,
                                cell.z_vector + cell.y_vector,
                                cell.x_vector + cell.z_vector,
                                cell.x_vector + cell.y_vector + cell.z_vector)
                        ],
                        *[
                            Cylinder(cell.base_vertex,
                                    cell.base_vertex + shift,
                                    material=cell_edge_mat,
                                    radius=0.2)
                            for shift in (cell.x_vector, cell.y_vector,
                                        cell.z_vector)
                        ],
                        *self.temp_intersections,
                    ]
        else:
                try:
                    result.append(
                        self.find_conflicts(
                            coords=p[index].xyz,
                            symbols=p.atom_symbols,
                            # For now, we test on fixed cell parameters
                            cell=self.start_cell_data,
                        ))
                except:
                    print(f'exception at frame {_index}, alpha {alpha}')
                    return None,_index,alpha
                save=False
                if save:
                    # r = MoleculeReconstructor(datafile=STARTING_DATAFILE, show_progress=True)
                    r = MoleculeReconstructor(datafile=self.datafile, show_progress=self.show_progress)
                    r.reconstruct_xyz(xyz_filename, xyz_filename[:-4]+f'-{index}'+'-check.xyz',p[index].xyz,p.atom_symbols)

        return result,_index,alpha

    def find_conflicts(
        self,
        coords: np.ndarray,
        symbols: list[str],
        cell: CellParameters,
    ) -> list[dict[str, tuple[int, ...]]]:
        """Find ring/bond intersections in a given conformation

        Args:
            coords (np.ndarray): Natoms x 3 matrix of atomic XYZ coordinatex
            symbols (list[str]): list of element symbols

        Returns:
            list[dict[str, tuple[int, ...]]]: list of intersection descriptions (see `analyze_xyz_file` for details)
        """
        self.crosscheck_symbols(symbols)

        if self.show_progress:
            ring_iter = tqdm.tqdm(self.ring_indices, leave=False)
        else:
            ring_iter = self.ring_indices
        crossing_bonds = []
        for i, cur_ring_indices in enumerate(ring_iter):
            # print(f"{i}/{len(ring_iter)}")
            # if i == 100:
            #     break
            # Verify that the ring is not ripped apart by some cell translation
            skip = False
            check_pairs = [(i, i + 1)
                           for i in range(len(cur_ring_indices) - 1)]
            check_pairs.append((0, len(cur_ring_indices) - 1))

            for atomA, atomB in check_pairs:
                bond_length_per = LAMMPSLoader.periodic_distance(
                    coords[cur_ring_indices[atomA]],
                    coords[cur_ring_indices[atomB]],
                    cell,
                )
                if bond_length_per > BOND_LENGTH_CUTOFF:
                    skip = True
                    break

            if skip:
                pass
                # raise Exception(f'Ring {cur_ring_indices} is ripped apart')

            crossing_bonds += self.find_crossing_bonds(
                coords=coords,
                ring_indices=cur_ring_indices,
                cell=cell,
            )
            # break

        return crossing_bonds

    def find_crossing_bonds(
        self,
        coords: np.ndarray,
        ring_indices: tuple[int, ...],
        cell: CellParameters,
    ) -> list[dict[str, tuple[int, ...]]]:
        """Identify ring/bond intersections for given atomic coordinates and ring atom indices.

        Args:
            coords (np.ndarray): XYZs of all atoms
            ring_indices (tuple[int, ...]): indices of atoms in some ring (indexing starts from 0)

        Returns:
            list[dict[str, tuple[int, ...]]]: list of intersection descriptions (see `analyze_xyz_file` for details)
        """
        assert self.symbols is not None
        assert self.molgraph is not None

        if self.debug_plotting:
            # A hack to pass coords to the 'fit_plane' method
            self.coords = coords

        ring_atoms_coords = coords[ring_indices, :]
        ring_frame = self.fit_plane(ring_atoms_coords, cell)

        # All radius vectors are computed relative to ring center
        coord_origin = ring_frame[:3, 3].copy()
        ring_frame[:3, 3] = np.zeros(3)

        sks_plane = Plane(point=np.zeros(3), normal=ring_frame[:3, 2])
        global_to_local = np.linalg.inv(ring_frame[:3, :3])

        adjusted_coords = np.array([
            LAMMPSLoader.mod_cell(xyz - coord_origin, cell) for xyz in coords
        ])
        ring_radius = max(
            [np.linalg.norm(v) for v in adjusted_coords[ring_indices, :]])

        node_accepted = {
            n: n not in ring_indices and np.linalg.norm(v) < ring_radius * 4.0
            for n, v in enumerate(adjusted_coords)
        }

        relative_plane_positions = [
            (global_to_local @ adjusted_coords[node])[2] > 0.0
            for node in range(self.molgraph.number_of_nodes())
        ]

        # Coordiantes of ring atom in plane of the ring
        ringatoms_local_coords = [
            (global_to_local @ adjusted_coords[atom_index])[:2]
            for atom_index in ring_indices
        ]

        # if self.debug_plotting:
        #     from visualization_pipelines.blender_scene import (
        #         SimplePoint,
        #         Cylinder,
        #         ArrowMaterial,
        #         CoordFrame,
        #     )
        #     assert self.symbols is not None

        #     cell_vertex_mat = ArrowMaterial(color='#FF0000',
        #                                     name='cell_vertex')
        #     cell_edge_mat = ArrowMaterial(color='#FFFF00', name='cell_edge')
        #     threed_mat = ArrowMaterial(
        #         color='#A65628',  # brown
        #         name='threed_mat')
        #     twod_mat = ArrowMaterial(
        #         color='#984EA3',  # purple
        #         name='twod_mat')
        #     centroid_mat = ArrowMaterial(
        #         color='#f781bf',  # pink
        #         name='centroid_mat')
        #     green_mat = ArrowMaterial(color='#00FF00', name='green_mat')
        #     blue_mat = ArrowMaterial(color='#0000FF', name='blue_mat')
        #     self.plot_data[f'intermid_{self.scene_id}'] = [
        #         # CELL
        #         *[
        #             SimplePoint(cell.base_vertex + shift,
        #                         material=cell_vertex_mat,
        #                         size=0.5) for shift in
        #             (np.zeros(3), cell.x_vector, cell.y_vector, cell.z_vector,
        #              cell.x_vector + cell.y_vector, cell.z_vector +
        #              cell.y_vector, cell.x_vector + cell.z_vector,
        #              cell.x_vector + cell.y_vector + cell.z_vector)
        #         ],
        #         *[
        #             Cylinder(cell.base_vertex,
        #                      cell.base_vertex + shift,
        #                      material=cell_edge_mat,
        #                      radius=0.2)
        #             for shift in (cell.x_vector, cell.y_vector, cell.z_vector)
        #         ],
        #         # DATA
        #         *[
        #             SimplePoint(
        #                 adjusted_coords[idx], material=threed_mat, size=0.3)
        #             for idx in ring_indices
        #         ],
        #         # *[
        #         #     SimplePoint(
        #         #         [*(xy.tolist()), 0.0], material=twod_mat, size=0.2)
        #         #     for xy in ringatoms_local_coords
        #         # ],
        #         *[
        #             SimplePoint(xyz,
        #                         material=green_mat
        #                         if relative_plane_positions[i] else blue_mat,
        #                         size=0.2)
        #             for i, xyz in enumerate(adjusted_coords)
        #             if np.linalg.norm(xyz) < 6.0
        #         ],
        #         # SimplePoint(coord_origin, material=centroid_mat, size=0.05),
        #         CoordFrame(ring_frame),
        #     ]
        # return []
        on_opposite_sides = lambda idxA, idxB: (relative_plane_positions[
            idxA] ^ relative_plane_positions[idxB])

        if self.debug_plotting:
            points = []
            plane_points = []

        detected_intersections = []
        ring_indices_plus_one = tuple(x + 1 for x in ring_indices)

        def save_intersection(A: int, B: int) -> None:
            """Store intersection info with indexing starting from 1
            """
            detected_intersections.append({
                'bond': (A + 1, B + 1),
                'ring': ring_indices_plus_one,
            })

        for nodeA, nodeB in self.molgraph.edges:
            if node_accepted[nodeA] and node_accepted[nodeB] and (
                    on_opposite_sides(nodeA, nodeB)):
                naive_direction = (adjusted_coords[nodeB] -
                                   adjusted_coords[nodeA])
                direction = LAMMPSLoader.mod_cell(naive_direction.copy(), cell)
                if (naive_direction != direction).any():
                    # Looks like these atoms are connected through cell walls and not the ring
                    continue
                if np.linalg.norm(direction) > BOND_LENGTH_CUTOFF:
                    pass
                    # raise Exception(
                    #     f"This bond is suspicious: {nodeA+1}-{nodeB+1}. "
                    #     f"Normalized length={np.linalg.norm(direction)}")

                sks_line = Line(point=adjusted_coords[nodeA],
                                direction=direction)
                sks_intersection_point = sks_plane.intersect_line(sks_line)
                intersection_point = np.array([*sks_intersection_point])
                intersection_point_projection = (
                    global_to_local @ intersection_point)[:2]

                if self.check_intersection_through_ring(
                        ring_atoms_xy=ringatoms_local_coords,
                        point_xy=intersection_point_projection,
                ):
                    if self.debug_plotting:
                        points.append((nodeA, nodeB, intersection_point,
                                       naive_direction, direction))
                        plane_points.append(
                            (nodeA, nodeB, intersection_point_projection))
                    save_intersection(nodeA, nodeB)

        # if self.debug_plotting and len(points) > 0:
        #     from visualization_pipelines.blender_scene import (
        #         SimplePoint,
        #         ArrowMaterial,
        #         Arrow,
        #     )
        #     assert self.symbols is not None

        #     threed_mat = ArrowMaterial(
        #         color='#A65628',  # brown
        #         name='threed_mat')
        #     twod_mat = ArrowMaterial(
        #         color='#984EA3',  # purple
        #         name='twod_mat')
        #     centroid_mat = ArrowMaterial(
        #         color='#f781bf',  # pink
        #         name='centroid_mat')
        #     green_mat = ArrowMaterial(color='#00FF00', name='green_mat')
        #     blue_mat = ArrowMaterial(color='#0000FF', name='blue_mat')
        #     self.temp_intersections += [
        #         *[
        #             SimplePoint(coords[idx], material=threed_mat, size=0.3)
        #             for idx in ring_indices
        #         ],
        #         *[
        #             point for nodeA, nodeB, intersection_xyz, naive_direction,
        #             direction in points for point in [
        #                 # SimplePoint(
        #                 #     intersection_xyz, material=centroid_mat, size=0.1),
        #                 SimplePoint(
        #                     coords[nodeA], material=green_mat, size=0.1),
        #                 SimplePoint(coords[nodeB], material=blue_mat,
        #                             size=0.1),
        #                 Arrow(coords[nodeA],
        #                       coords[nodeA] + direction,
        #                       material=twod_mat),
        #             ]
        #         ],
        #     ]
        #     # self.plot_data[
        #     #     f'intersect-{self.scene_id}-{"*".join(str(i) for i in ring_indices)}'] = [
        #     #         *[
        #     #             SimplePoint(cell.base_vertex + shift,
        #     #                         material=cell_vertex_mat,
        #     #                         size=0.5)
        #     #             for shift in (
        #     #                 np.zeros(3), cell.x_vector, cell.y_vector,
        #     #                 cell.z_vector, cell.x_vector + cell.y_vector,
        #     #                 cell.z_vector + cell.y_vector,
        #     #                 cell.x_vector + cell.z_vector,
        #     #                 cell.x_vector + cell.y_vector + cell.z_vector)
        #     #         ],
        #     #         *[
        #     #             Cylinder(cell.base_vertex,
        #     #                      cell.base_vertex + shift,
        #     #                      material=cell_edge_mat,
        #     #                      radius=0.2)
        #     #             for shift in (cell.x_vector, cell.y_vector,
        #     #                           cell.z_vector)
        #     #         ],
        #     #         *[
        #     #             SimplePoint(adjusted_coords[idx],
        #     #                         material=threed_mat,
        #     #                         size=0.3) for idx in ring_indices
        #     #         ],
        #     #         *[
        #     #             point for nodeA, nodeB, intersection_xyz,
        #     #             naive_direction, direction in points for point in [
        #     #                 SimplePoint(intersection_xyz,
        #     #                             material=centroid_mat,
        #     #                             size=0.1),
        #     #                 SimplePoint(adjusted_coords[nodeA],
        #     #                             material=green_mat,
        #     #                             size=0.1),
        #     #                 SimplePoint(adjusted_coords[nodeB],
        #     #                             material=blue_mat,
        #     #                             size=0.1),
        #     #                 Arrow(adjusted_coords[nodeA],
        #     #                       adjusted_coords[nodeA] + naive_direction,
        #     #                       material=green_mat),
        #     #                 Arrow(adjusted_coords[nodeA],
        #     #                       adjusted_coords[nodeA] + direction,
        #     #                       material=twod_mat),
        #     #             ]
        #     #         ],
        #     #         # CoordFrame(ring_frame),
        #     #     ]

        return detected_intersections

    def fit_plane(
        self,
        points_raw: np.ndarray,
        cell: CellParameters,
    ) -> np.ndarray:
        """Identifies the plane of ring. Result is returned as coordinate frame
        at the center of ring and Z axis perpendicular to the ring.
        Coordinate frame is represented via SE(3) matrix (see https://doi.org/10.1002/wcms.1690)

        Args:
            points_raw (np.ndarray): Coordinates of atoms in ring (in global coordinate frame).

        Returns:
            np.ndarray: 4x4 SE(3) matrix of coordinate frame in the plane of ring
        """

        # Choose the first point as a coordinate origin
        origin = points_raw[0]
        points = [
            LAMMPSLoader.mod_cell(point - origin, cell) for point in points_raw
        ]

        # Center the points by subtracting the centroid
        centroid = np.mean(points, axis=0)
        centered_points = points - centroid
        _, _, Vt = np.linalg.svd(centered_points)

        frame = np.zeros((4, 4))
        frame[:3, :3] = Vt.T
        frame[:3, 3] = centroid + origin
        frame[3, 3] = 1

        # if self.debug_plotting:
        #     from visualization_pipelines.blender_scene import (
        #         SimplePoint,
        #         Cylinder,
        #         ArrowMaterial,
        #         CoordFrame,
        #     )
        #     assert self.symbols is not None

        #     cell_vertex_mat = ArrowMaterial(color='#FF0000',
        #                                     name='cell_vertex')
        #     cell_edge_mat = ArrowMaterial(color='#FFFF00', name='cell_edge')
        #     raw_point_mat = ArrowMaterial(color='#A65628',
        #                                   name='raw_point_mat')
        #     centroid_mat = ArrowMaterial(color='#f781bf', name='centroid_mat')
        #     point_mat = ArrowMaterial(color='#984EA3', name='point_mat')
        #     sc = 1.0
        #     self.plot_data[f'plane_{self.scene_id}'] = [
        #         *[
        #             SimplePoint(sc * (cell.base_vertex + shift),
        #                         material=cell_vertex_mat,
        #                         size=0.5) for shift in
        #             (np.zeros(3), cell.x_vector, cell.y_vector, cell.z_vector,
        #              cell.x_vector + cell.y_vector, cell.z_vector +
        #              cell.y_vector, cell.x_vector + cell.z_vector,
        #              cell.x_vector + cell.y_vector + cell.z_vector)
        #         ],
        #         *[
        #             Cylinder(sc * cell.base_vertex,
        #                      sc * (cell.base_vertex + shift),
        #                      material=cell_edge_mat,
        #                      radius=0.2)
        #             for shift in (cell.x_vector, cell.y_vector, cell.z_vector)
        #         ],
        #         *[
        #             SimplePoint(sc * xyz, material=raw_point_mat, size=0.05)
        #             for xyz in points_raw
        #         ],
        #         *[
        #             SimplePoint(sc * xyz, material=point_mat, size=0.05)
        #             for xyz in points
        #         ],
        #         SimplePoint(sc * centroid, material=centroid_mat, size=0.05),
        #         CoordFrame(frame),
        #     ]
        return frame

    def check_intersection_through_ring(
        self,
        point_xy: np.ndarray,
        ring_atoms_xy: list[np.ndarray],
    ) -> bool:
        """Check whether the point `point_xy` is inside a polygon formed by points `ring_atoms_xy`.
        It verifies that the intersection point of the line connecting two atoms is inside the ring.

        Args:
            point_xy (np.ndarray): 2D coordinates of intersection point in the plane of ring
            ring_atoms_xy (list[np.ndarray]): Coordinates of atoms in the plane of ring

        Returns:
            bool: True if bond DOES go through the ring, False if it does not
        """
        from scipy.spatial import ConvexHull
        hull = ConvexHull([point_xy] + ring_atoms_xy)
        return 0 not in hull.vertices


class MoleculeReconstructor(LAMMPSLoader):
    SOURCE_NODE = 0  # The choice of the source does not matter. Nodes are indexed from 0 onwards

    def graph_preprocessing(self) -> None:
        assert self.molgraph is not None
        print('Doing topology analysis...')
        self.assembly_graph, self.node_order, self.conn_components = self.build_assembly_graph(
            self.molgraph)
        print('Finished topology analysis')

    def build_assembly_graph(
            self, input_graph: nx.Graph
    ) -> tuple[nx.Graph, list[int], list[list[int]]]:
        assert self.molgraph is not None
        assert input_graph.number_of_nodes() != 0

        order_of = {}
        node_order_full = []
        components = [list(c) for c in nx.connected_components(input_graph)]
        for component_nodes in components:
            molsubgraph = input_graph.subgraph(component_nodes)
            node_order = list(
                nx.dfs_preorder_nodes(molsubgraph,
                                      source=next(iter(component_nodes))))
            order_of.update({
                node: order
                for order, node in enumerate(node_order, start=len(order_of))
            })
            node_order_full += node_order

        res_graph = nx.DiGraph()
        res_graph.add_nodes_from(input_graph.nodes)
        for nA, nB in input_graph.edges:
            if order_of[nA] < order_of[nB]:
                res_graph.add_edge(nA, nB)
            else:
                res_graph.add_edge(nB, nA)

        return res_graph.reverse(), node_order_full, components

    def reconstruct_xyz(self, input_path: str, output_path: str, xyz=None,sym=None) -> None:
        """Reconstruct molecules in XYZ file that are messed up due to unit cell translations

        Args:
            input_path (str): input XYZ-file where all atoms are forced inside unit cell
            output_path (str): output XYZ-file of intelligible geometries reconstructed based on input connectivity
        """
        input_p = Confpool()
        if xyz is None:
            input_p.include_from_file(input_path)
        else:
            input_p.include_from_xyz(xyz,'--')
            input_p.atom_symbols = sym
            
        self.crosscheck_symbols(input_p.atom_symbols)

        output_p = Confpool()
        if self.show_progress:
            pool_iter = enumerate(tqdm.tqdm(input_p))
        else:
            pool_iter = enumerate(input_p)
        for i, m in pool_iter:
            fixed_xyz: np.ndarray = self.reconstruct_single_structure(
                m.xyz.copy(),
                # For now, we test on fixed cell parameters
                cell=self.start_cell_data,
            )
            if fixed_xyz is None:
                assert self.debug_plotting
                return None

            output_p.include_from_xyz(fixed_xyz, m.descr)

        output_p.atom_symbols = input_p.atom_symbols
        # output_p.save_xyz(output_path)
        output_p.save(output_path)

    def reconstruct_single_structure(
        self,
        coords: np.ndarray,
        cell: CellParameters,
    ) -> np.ndarray:
        """WARNING: the `coords` array gets modified here
        """
        # if self.debug_plotting:

        #     def dump_scene(i, prev_coords, new):
        #         from visualization_pipelines.blender_scene import (
        #             SimplePoint,
        #             Cylinder,
        #             ArrowMaterial,
        #             CoordFrame,
        #         )
        #         ic(i, prev_coords)
        #         # raise Exception("AAAAAAAAAAA")
        #         assert self.symbols is not None

        #         cell_vertex_mat = ArrowMaterial(color='#FF0000',
        #                                         name='cell_vertex')
        #         cell_edge_mat = ArrowMaterial(color='#FFFF00',
        #                                       name='cell_edge')
        #         raw_point_mat = ArrowMaterial(color='#A65628',
        #                                       name='raw_point_mat')
        #         centroid_mat = ArrowMaterial(color='#f781bf',
        #                                      name='centroid_mat')
        #         point_mat = ArrowMaterial(color='#984EA3', name='point_mat')
        #         green_mat = ArrowMaterial(color='#00FF00', name='green_mat')
        #         blue_mat = ArrowMaterial(color='#0000FF', name='blue_mat')
        #         self.plot_data[f'reconstr'] = [
        #             *[
        #                 SimplePoint(cell.base_vertex + shift,
        #                             material=cell_vertex_mat,
        #                             size=0.5)
        #                 for shift in (
        #                     np.zeros(3), cell.x_vector, cell.y_vector,
        #                     cell.z_vector, cell.x_vector + cell.y_vector,
        #                     cell.z_vector + cell.y_vector,
        #                     cell.x_vector + cell.z_vector,
        #                     cell.x_vector + cell.y_vector + cell.z_vector)
        #             ],
        #             *[
        #                 Cylinder(cell.base_vertex,
        #                          cell.base_vertex + shift,
        #                          material=cell_edge_mat,
        #                          radius=0.2)
        #                 for shift in (cell.x_vector, cell.y_vector,
        #                               cell.z_vector)
        #             ],
        #             *[
        #                 SimplePoint(
        #                     coords[node], material=raw_point_mat, size=0.2)
        #                 for node in self.node_order[:i]
        #             ],
        #             *[
        #                 SimplePoint(xyz, material=green_mat, size=0.4)
        #                 for xyz in prev_coords
        #             ],
        #             SimplePoint(new, material=blue_mat, size=0.4),
        #             *[
        #                 Cylinder(coords[nA],
        #                          coords[nB],
        #                          material=point_mat,
        #                          radius=0.1)
        #                 for nA, nB in self.molgraph.edges if
        #                 nA in self.node_order[:i] and nB in self.node_order[:i]
        #             ],
        #         ]

        for i, cur_node in enumerate(self.node_order):
            new_coords = []
            for parent_node in self.assembly_graph.neighbors(cur_node):
                new = coords[parent_node] + self.mod_cell(
                    coords[cur_node] - coords[parent_node], cell)

                assert all(
                    np.linalg.norm(new - prev) < 0.001 for prev in new_coords
                ), (f"Unable to reconstruct some cycle involving bond {cur_node}-{parent_node}. "
                    f"Looks like this is impossible: {new} vs. {repr(new_coords)}"
                    )
                # if self.debug_plotting and not all(
                #         np.linalg.norm(new - prev) < 0.001
                #         for prev in new_coords):
                #     dump_scene(i, new_coords, new)
                #     return None
                new_coords.append(new)

            if len(new_coords) > 0:
                coords[cur_node] = new_coords[0]

        for component in self.conn_components:
            comp_xyz = coords[component, :]
            centroid_init = comp_xyz.mean(axis=0)
            # ic(max(comp_xyz, key=lambda v: np.linalg.norm(v)))
            shift = self.mod_cell(centroid_init, cell) - centroid_init
            comp_xyz += shift

        return coords


if __name__ == "__main__":
    d = CrossingDetector(datafile=STARTING_DATAFILE, show_progress=True)
    # with cProfile.Profile() as pr:
    result, *_ = d.analyze_xyz_file(STARTING_XYZ,index=INDEX)
    # pr.dump_stats('numba.pstats')
    # print(result)
    for i, overlap_data in enumerate(result, start=1):
        # overlap_data, *_ = res
        if len(overlap_data) == 0:
            print(f"{i}) No ring crossings")
        else:
            print(f"{i}) Ring crossings found: {overlap_data}")
            res = []
            for d in overlap_data:
                for v in d.values():
                    for elem in v:
                        res.append(elem)
            print(f'Total number of crossings: {len(overlap_data)}')
            print(f'Atom selection for VMD: {" ".join([str(i) for i in set(res)])}')

    # r = MoleculeReconstructor(datafile=STARTING_DATAFILE, show_progress=True)
    # if index !=None:
        
    #     r.reconstruct_xyz(STARTING_XYZ, STARTING_XYZ[:-4]+'-check.xyz')
    # else:
    #     r.reconstruct_xyz(STARTING_XYZ, STARTING_XYZ[:-4]+'-check.xyz')