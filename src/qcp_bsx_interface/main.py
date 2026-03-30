# from augmentation import requires_augmentation, produce_augmentation
# If basis does not require augmentation, just use the original basis set

# from util.import_helper import *
# 
# recognized_suffixes=['wc3','c3'] # make core polarization for third row elements, or all elements respectively
# 
# def has_custom_modification(basis_name):
#     parts=basis_name.split('-')
#     if parts[-1].lower() in [ x.lower() for x in recognized_suffixes]:
#         return True
#     return False
# 
# def requires_modification(basis_name):
#     return has_custom_modification(basis_name) or requires_augmentation(basis_name)
# 
# def create_basis(basis_name, elements):
#     """ The basis may be constituted by multiple segments
#     pref-(aug-cc-pvqz)-suffix
#     d-aug, t-aug for augmentation
#     """
#     aug_req=requires_augmentation(basis_name)
#     custom_modification=has_custom_modification(basis_name)
# 
#     prefix=basis_name.split('-aug-')[0]
#     assert prefix in ['d','t','q']
# 
#     if aug_req:
#         elements = produce_augmentation(elements)
#     return elements