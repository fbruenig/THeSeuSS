'''Preparation and modification of input files'''

#AUTHOR: Ariadni Boziki

import re


class GeometryProcessor():

    """
    Class for exctracting information out of a geometry input file.

    Attributes
    ----------

    geom_input: str
        Namefile of the geometry input file.
    code: str
        Code
    """

    def __init__(self, geom_input: str, code: str):

        self.geom_input = geom_input
        self.code = code

    def number_of_atoms(self)-> int:
        """
        It returns the number of atoms in the structure.

        index: int
            Number of atoms in a given geometry input file.
        """

        index = 0
        with open(self.geom_input,'r') as geometry_input:
            for lines in geometry_input:
                if self.code == 'aims' and (("atom" in lines or "atom_frac" in lines) and "hessian_block_lv_atom" not in lines):
                    index += 1
                elif self.code == 'dftb+':
                    no_of_columns = len(lines.split())
                    if no_of_columns == 5:
                        index += 1
        return index

    def read_lattice(self)-> list:
        """
        It returns the lattice cell (9 elements) from the geometry input file.

        cell: list
            List containing the lattice cell elements as floats.
        """
	
        cell = []

        with open(self.geom_input, 'r') as geometry_input:

            if self.code == 'aims':
                for lines in geometry_input:
                    if 'lattice_vector' in lines:
                        cell_tmp = lines.split()[1:]
                        cell_tmp_float = [float(s.replace(',','')) for s in cell_tmp]
                        cell.append(cell_tmp_float)

            elif self.code == 'dftb+':
                for line in (geometry_input.readlines() [-3:]):
                    cell_tmp = line.split()
                    cell_tmp_float = [float(s.replace(',','')) for s in cell_tmp]
                    cell.append(cell_tmp_float)
        return cell

    def read_coordinates(self)-> list:
        """
        It returns the coordinates from the geometry input file.
        
        coord : list
            List containing coordinate tuples as floats.
        """

        coord = []	

        with open(self.geom_input, 'r') as geometry_input:

            if self.code == 'aims':
                for lines in geometry_input:
                    if ("atom" in lines or "atom_frac" in lines) and "hessian_block_lv_atom" not in lines:
                        coord_tmp = lines.split()[1:4]
                        coord_tmp_float = [float(s.replace(',','')) for s in coord_tmp]
                        coord.append(coord_tmp_float)

            elif self.code == 'dftb+':
                for lines in geometry_input:
                    no_of_columns = len(lines.split())
                    if no_of_columns == 5:
                        coord_tmp = lines.split()[2:5]
                        coord_tmp_float = [float(s.replace(',','')) for s in coord_tmp]
                        coord.append(coord_tmp_float)
        return coord

    def read_the_atom_type(self)-> [list, list]:
        """
        It returns the atom type of the atoms of the structure from the geometry input file.

        atom_type_no: tuple
            Tuple containing lists of atom type numbers.
        atom_type: tuple
            Tuple containing lists of atom types.
        """

        atom_type = []
        atom_type_no = []
        atom_types_tmp_dftb = []

        with open(self.geom_input,'r') as geometry_input:

            contents = geometry_input.readlines()

            if self.code == 'dftb+':
                atom_types_tmp_dftb = contents[1].split()
                
            for lines in contents:
                
                if self.code == 'aims' and (("atom" in lines or "atom_frac" in lines) and "hessian_block_lv_atom" not in lines):
                    atom_type_tmp = lines.split()[-1]
                    atom_type.append(atom_type_tmp)
                    
                    if len(atom_type) == 1:
                        atom_type_no_tmp = 1
                    elif atom_type[-1] == atom_type[-2]:
                        atom_type_no_tmp = atom_type_no[-1]
                    else:
                        atom_type_no_tmp = atom_type_no[-1] + 1
                    
                    atom_type_no.append(atom_type_no_tmp)
                
                elif self.code == 'dftb+':
                    no_of_columns = len(lines.split())
                    if no_of_columns == 5:
                        atom_type_no_tmp = lines.split()[1]
                        atom_type_no.append(atom_type_no_tmp)
            if (self.code == 'dftb+'):
                for ii in range(len(atom_types_tmp_dftb)):
                    kk = ii+1
                    for jj in atom_type_no:
                        ll = int(jj)
                        if kk is ll:
                            atom_type.append(atom_types_tmp_dftb[ii])
        return atom_type_no, atom_type
