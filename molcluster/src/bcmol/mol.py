from __future__ import annotations

import copy
from pathlib import Path
from typing import Iterator, Literal, TextIO

from .constants import ATOMIC_NUMBERS, ELEMENTS
from .core import Atom, Bond, BondType, Chain
from .core import Molecule as BaseMolecule
from .core import Residue, ResType, guess_mol_type, guess_residue_type


def parse_bonds(mol: Molecule, bonds: list[tuple[int, int, BondType]]):
    for bond in bonds:
        atom, partner = mol.atoms[bond[0] - 1], mol.atoms[bond[1] - 1]
        mol.bonds.append(Bond(atom=atom, partner=partner, type=bond[2]))


class Molecule(BaseMolecule):
    def clone(self) -> Molecule:
        return copy.deepcopy(self)

    def select(self, atoms: list[int]) -> Molecule:
        mol = self.clone()

        mol.atoms = [x for x in mol.atoms if x.idx in atoms]
        mol.bonds = [x for x in mol.bonds if
                     x.atom in mol.atoms and x.partner in mol.atoms]

        residues = set(x.residue.idx for x in mol.atoms)
        mol.residues = [x for x in mol.residues if x.idx in residues]

        chains = set(x.chain.name for x in mol.residues)
        mol.chains = [x for x in mol.chains if x.name in chains]

        for i, atom in enumerate(mol.atoms):
            atom.idx = i
        for i, residue in enumerate(mol.residues):
            residue.idx = i
        for i, chain in enumerate(mol.chains):
            chain.idx = i

        return mol

    def write(self,
              file: str | Path = None,
              format: Literal['mol2', 'sdf', 'pdb'] = None,
              idx: int = None,
              increase_tercount: bool = True,
              add_end_card: bool = True,
              write_conects: bool = True,
              ) -> str | None:
        """

        Args:
            file (str | Path): Output file. If None, return the content instead
                    of writing to file.
            format (str): Output molecular format.
            idx (int): Model index which starts from 1.
            increase_tercount (bool): Whether to increase the atom number for
                    PDB ter cards.
            add_end_card (bool): Whether to add PDB `END` card. Sometimes it is
                    necessary not to add the end card for intermediate models.
            write_conects (bool): Whether to write PDB connectivity records.

        Returns:

        """
        if not file and not format:
            raise Exception('请提供输出文件名或输出格式！')

        if file:
            if isinstance(file, str):
                file = Path(file)
            format = format or file.suffix[1:]

        from .mol2 import Writer as Mol2Writer
        from .pdb import Writer as PDBWriter
        from .sdf import Writer as SDFWriter

        match format:
            case 'mol2':
                writer = Mol2Writer()
            case 'pdb':
                writer = PDBWriter(increase_tercount)
            case 'sdf':
                writer = SDFWriter()
            case _:
                raise Exception(f"Unsupported format: {format}!")

        if format == 'pdb':
            content = writer.write(self, idx, write_conects=write_conects)
            if add_end_card:
                content += 'END\n'
        else:
            content = writer.write(self)

        if file:
            file.write_text(content)
        else:
            return content

    def to_pybel(self):
        from openbabel import pybel

        # Create openbabel Molecule object
        m = pybel.ob.OBMol()

        # Add residues & chains
        residues = []
        for residue in self.residues:
            r = m.NewResidue()
            r.SetName(residue.name)
            r.SetNum(residue.number)
            r.SetChain(residue.chain.name)
            residues.append(r)

        # Add atoms
        for atom in self.atoms:
            a = m.NewAtom()
            a.SetAtomicNum(ATOMIC_NUMBERS[atom.element])
            a.SetVector(atom.x, atom.y, atom.z)
            r = next(r for r in residues
                     if r.GetNum() == atom.residue.number
                     and r.GetName() == atom.residue.name
                     and r.GetChain() == atom.residue.chain.name)
            a.SetResidue(r)
            if atom.charge is not None:
                a.SetPartialCharge(atom.charge)
            a.SetTitle(atom.name)
            a.SetType(atom.type)

        # Add bonds
        for bond in self.bonds:
            m.AddBond(bond.atom.idx + 1, bond.partner.idx + 1, bond.order)

        m = pybel.Molecule(m)
        m.title = self.name
        for key, value in self.data.items():
            m.data[key] = value

        return m


def from_pybel(mol) -> Molecule:
    from openbabel import pybel

    m = Molecule(name=mol.title,
                 comment='',
                 charge=mol.charge,
                 charge_type=dict(mol.data).get('PartialCharges', 'NO_CHARGES'),
                 data={k: v for k, v in mol.data.items()
                       if k not in ('PartialCharges', 'OpenBabel Symmetry Classes')},
                 )

    # Make residues & chains
    chains_dict = {}
    for residue in mol.residues:
        chain_name = residue.OBResidue.GetChain()
        if chain := chains_dict.get(chain_name) is None:
            chains_dict[chain_name] = Chain(idx=0, name=chain_name)

        res = Residue(idx=residue.idx,
                      name=residue.name[:3],
                      number=residue.OBResidue.GetNum(),
                      type=guess_residue_type(residue.name[:3]),
                      chain=chain)
        m.residues.append(res)

    for i, chain in enumerate(sorted(chains_dict.values(), key=lambda x: x.name)):
        chain.idx = i
        m.chains.append(chain)

    # Guess mol type
    m.type = guess_mol_type(m)

    # Prepare for atom types translation
    # https://openbabel.org/dev-api/classOpenBabel_1_1OBTypeTable.shtml
    ttab = pybel.ob.OBTypeTable()
    ttab.SetFromType('INT')
    ttab.SetToType('SYB')

    # Make atoms
    residues_dict = {residue.idx: residue for residue in m.residues}
    for i, atom in enumerate(mol.atoms):
        residue = residues_dict[atom.residue.OBResidue.GetIdx()]
        element = ELEMENTS.get(atom.atomicnum, 'X')
        a = Atom(idx=i,
                 number=atom.idx,
                 name=atom.residue.OBResidue.GetAtomID(atom.OBAtom),
                 # https://openbabel-discuss.narkive.com/4A835GV1/open-babel-get-atom-name-information-from-obmol-class-using-python-bindings
                 x=atom.coords[0],
                 y=atom.coords[1],
                 z=atom.coords[2],
                 element=element,
                 type=ttab.Translate(atom.OBAtom.GetType()),
                 charge=atom.partialcharge,
                 is_het=residue.type in (ResType.hetero, ResType.small, ResType.unknown),
                 residue=residue,
                 )
        m.atoms.append(a)

    # Make bonds
    bonds = []
    for i in range(mol.OBMol.NumBonds()):
        bond = mol.OBMol.GetBondById(i)
        begin_idx, end_idx = bond.GetBeginAtom().GetIdx(), bond.GetEndAtom().GetIdx()
        if begin_idx > end_idx:
            begin_idx, end_idx = end_idx, begin_idx
        bonds.append((begin_idx, end_idx, BondType.get_bond_type_by_count(bond.GetBondOrder())))

    bonds.sort()
    parse_bonds(m, bonds)

    return m


def to_pybel(mol: Molecule):
    from openbabel import pybel

    # Create openbabel Molecule object
    m = pybel.ob.OBMol()

    # Add residues & chains
    residues = []
    for residue in mol.residues:
        r = m.NewResidue()
        r.SetName(residue.name)
        r.SetNum(residue.number)
        r.SetChain(residue.chain.name)
        residues.append(r)

    # Add atoms
    for atom in mol.atoms:
        a = m.NewAtom()
        a.SetAtomicNum(ATOMIC_NUMBERS[atom.element])
        a.SetVector(atom.x, atom.y, atom.z)
        r = next(r for r in residues
                 if r.GetNum() == atom.residue.number
                 and r.GetName() == atom.residue.name
                 and r.GetChain() == atom.residue.chain.name)
        a.SetResidue(r)
        a.SetPartialCharge(atom.charge)
        a.SetTitle(atom.name)
        a.SetType(atom.type)

    # Add bonds
    for bond in mol.bonds:
        m.AddBond(bond.atom.idx + 1, bond.partner.idx + 1, bond.order)

    m = pybel.Molecule(m)
    m.title = mol.name
    for key, value in mol.data.items():
        m.data[key] = value

    return m


def read_file(file: str | Path | TextIO, format: str = None) -> Iterator[Molecule]:
    from .mol2 import Reader as Mol2Reader
    from .pdb import Reader as PDBReader
    from .sdf import Reader as SDFReader

    if isinstance(file, str):
        file = Path(file)
    if isinstance(file, Path):
        fio = file.open()
        name = file.stem
        format = file.suffix[1:]
    else:
        fio = file
        name = 'MOL'

    match format:
        case 'mol2':
            reader = Mol2Reader(name=name)
        case 'pdb':
            reader = PDBReader(name=name)
        case 'sdf':
            reader = SDFReader(name=name)
        case _:
            raise Exception(f"Unsupported format: {format}!")

    return reader.read(fio)


def write_file(mols: Molecule | list[Molecule],
               file: str | Path = None,
               format: str = None,
               increase_tercount: bool = True) -> str | None:
    if not file and not format:
        raise Exception('请提供输出文件名或输出格式！')

    if file:
        if isinstance(file, str):
            file = Path(file)
        format = format or file.suffix[1:]

    if isinstance(mols, Molecule):
        mols = [mols]

    if len(mols) == 1:
        content = mols[0].write(format=format,
                                increase_tercount=increase_tercount,
                                add_end_card=True)
    else:
        content = ''
        for i, mol in enumerate(mols, 1):
            content += mol.write(format=format,
                                 idx=i,
                                 increase_tercount=increase_tercount,
                                 add_end_card=False,
                                 write_conects=i == len(mols),
                                 )
        if format == 'pdb':
            content += 'END\n'

    if file:
        file.write_text(content)
    else:
        return content


class OutputFile:
    def __init__(self,
                 file: str | Path | TextIO,
                 format: str = None,
                 increase_tercount: bool = True):
        if isinstance(file, str):
            file = Path(file)
        if isinstance(file, Path):
            format = format or file.suffix[1:]
            file = file.open('w')

        self.format = format
        self.file = file
        self.increase_tercount = increase_tercount
        self.last_mol = None

    def write(self, mol: Molecule, idx: int = None):
        content = mol.write(self.format,
                            idx=idx,
                            increase_tercount=self.increase_tercount,
                            add_end_card=False)
        self.file.write(content)
        self.last_mol = mol

    def close(self):
        if self.format == 'pdb':
            from .pdb import Writer

            writer = Writer(self.increase_tercount)
            self.file.write(writer.connect_records(self.last_mol.bonds))
            self.file.write('END\n')

        self.file.close()
