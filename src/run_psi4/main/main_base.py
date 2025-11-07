# def xyz_file_read(file_name):
#     with open(file_name, 'r') as rd:
#         lines=rd.readlines()
#     assert len(lines)>2
#     coordinates=[ x.split()[1:4] for x in lines[2:] if x.strip()!='' ]
#     atom_types=[ x.split()[0] for x in lines[2:] if x.strip()!='' ]
#     num_atoms=int(lines[0])
#     assert len(coordinates)==num_atoms, f"{len(coordinates)} {num_atoms}"
#     coordinates=np.array([ list(map(float, x)) for x in coordinates])
# 
#     return atom_types, coordinates
# 
# def xyz_to_string(file_name):
#     """
#     """
#     atom_types, coordinates = xyz_file_read(file_name)
#     molecule_str=molecule_string(atom_types=atom_types, coordinates=coordinates, units={'LENGTH':'ANGSTROM'})
# 
#     return molecule_str






