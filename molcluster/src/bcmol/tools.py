from dataclasses import fields
from types import GenericAlias

from .mol import *


def merge_mols(mols: list[Molecule], continue_residue_numbers: bool = True) -> Molecule:
    """Merge multiple molecules into a single one. Atoms are auto renumbered.
     Residue numbers will be kept if `continue_residue_numbers` is False; Otherwise,
     subsequent residues will be renumbered continuously.
     """
    if not mols:
        raise Exception('No molecules provided!')
    elif len(mols) == 1:
        return mols[0].clone()

    com = Molecule(name='Complex', comment='')

    chain_names = set()
    last_residue_number = 0
    for mol in mols:
        mol = mol.clone()

        # Make chains
        for chain in mol.chains:
            if chain.name not in chain_names:
                chain_names.add(chain.name)
                com.chains.append(Chain(idx=0, name=chain.name, mol=com))

        # Make residues
        mol.residues[-1].ter = True
        for i, residue in enumerate(mol.residues, last_residue_number + 1):
            chain = next(x for x in com.chains if x.name == residue.chain.name)
            residue.chain = chain
            if continue_residue_numbers:
                residue.number = i

        com.residues.extend(mol.residues)
        last_residue_number = com.residues[-1].number

        # Make atoms
        com.atoms.extend(mol.atoms)

        # Make bonds
        com.bonds.extend(mol.bonds)

        # Set data
        com.data |= mol.data

    # Reindex & renumber
    for i, residue in enumerate(com.residues):
        residue.idx = i

    for i, atom in enumerate(com.atoms):
        atom.idx = i
        atom.number = i + 1
    com.bonds.sort(key=lambda x: (x.atom.idx, x.partner.idx))

    # Set charge type & charge
    charge_types = {x.charge_type for x in mols}
    charges = [x.charge for x in mols]
    if 'NO_CHARGES' in charge_types or None in charges:
        com.charge_type = 'NO_CHARGES'
        com.charge = None
    else:
        com.charge = int(sum(charges))
        if len(charge_types) >= 2:
            com.charge_type = 'USER_CHARGES'
        else:
            com.charge_type = charge_types.pop()

    # Guess mol type
    com.type = guess_mol_type(com)

    return com


def check_bonds(mol: Molecule):
    bonds = [(x.atom.number, x.partner.number) for x in mol.bonds]
    for atom in mol.atoms:
        if atom.partners:
            for partner in atom.partners:
                pair = (atom.number, partner.number)
                if pair not in bonds and (partner.number, atom.number) not in bonds:
                    print(pair)


def remove_duplicate_objects(objects: list) -> list:
    if not objects:
        return []

    fields_ = tuple(f.name for f in fields(objects[0]) if not isinstance(f.type, GenericAlias)
                    and f.type not in ('Chain', 'Residue', 'Atom', 'Bond'))
    items = [[x, tuple(getattr(x, f) for f in fields_)] for x in objects]

    objects_ = []
    seen = set()
    for x, item in items:
        if item not in seen:
            seen.add(item)
            objects_.append(x)
    return objects_
