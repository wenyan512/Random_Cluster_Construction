from __future__ import annotations

import copy
from dataclasses import dataclass, field, replace
from enum import Enum

from .constants import RESTYPES


def deepcopy(x, memo=None):
    if memo is None:
        memo = {}

    y = memo.get(id(x))
    if y is not None:
        return y
    return copy.deepcopy(x, memo)


class BondType(str, Enum):
    zero = 'nc'
    single = '1'
    double = '2'
    triple = '3'
    amide = 'am'
    aromatic = 'ar'
    dummy = 'du'
    dative = 'da'
    unknown = 'un'

    @staticmethod
    def get_bond_type_by_count(count: int):
        bond_type_map = {
            1: BondType.single,
            2: BondType.double,
            3: BondType.triple,
        }
        return bond_type_map.get(count, BondType.unknown)


@dataclass
class Bond:
    atom: Atom = None
    partner: Atom = None
    type: BondType = BondType.single
    order: int = 1

    def __post_init__(self):
        order_map = {
            BondType.zero: 0,
            BondType.single: 1,
            BondType.double: 2,
            BondType.triple: 3,
            BondType.amide: 1,
            BondType.aromatic: 1,
            BondType.dummy: 1,
            BondType.dative: 1,
        }
        self.order = order_map.get(self.type, 1)

    def __repr__(self):
        return (f"Bond(atom=<Atom {self.atom.number}.{self.atom.name}>, "
                f"partner=<Atom {self.partner.number}.{self.partner.name}>, "
                f"type='{self.type}', order={self.order})")


@dataclass
class Atom:
    idx: int
    number: int
    name: str
    x: float
    y: float
    z: float
    element: str
    type: str = ''
    alt: str = ''
    insertion: str = ''
    segid: str = ''
    occupancy: float = 1.
    bfactor: float = 0.
    charge: float = None
    is_het: bool = False
    residue: Residue | None = None

    def __repr__(self):
        partners = ', '.join([f'<Atom {x.number}.{x.name}>' for x in self.partners])
        return (f"Atom(idx={self.idx}, number={self.number}, name='{self.name}', "
                f"x={self.x}, y={self.y}, z={self.z}, element='{self.element}', type='{self.type}', "
                f"alt='{self.alt}', insertion='{self.insertion}', segid='{self.segid}', "
                f"occupancy={self.occupancy}, bfactor={self.bfactor}, charge={self.charge}, "
                f"is_het={self.is_het}, residues=<Residue {self.residue.name}{self.residue.number}>, partners=[{partners}]")

    def __deepcopy__(self, memo=None):
        if memo is None:
            memo = {}

        if atom := memo.get(id(self)):
            return atom

        atom = replace(self, residue=None)
        memo[id(self)] = atom
        atom.residue = deepcopy(self.residue, memo)

        return atom

    @property
    def partners(self) -> list[Atom]:
        atoms = []
        for bond in self.residue.chain.mol.bonds:
            if bond.atom == self:
                atoms.append(bond.partner)
            elif bond.partner == self:
                atoms.append(bond.atom)
        return atoms


class ResType(str, Enum):
    bio = 'Biopolymer'
    protein = 'Protein'
    dna = 'DNA'
    rna = 'RNA'
    sugar = 'Saccharide'
    hetero = 'Hetero'
    small = 'Small'
    metal = 'Metal'
    ion = 'Ion'
    water = 'Water'
    unknown = 'Unknown'


@dataclass
class Residue:
    idx: int
    name: str
    number: int
    type: ResType = None
    chain: Chain | None = None
    ter: bool = False
    atoms_: list[Atom] = field(default_factory=list)

    def __repr__(self):
        atoms = ', '.join([f'<Atom {x.number}.{x.name}>' for x in self.atoms])
        return (f"Residue(name='{self.name}', number={self.number}, type='{self.type}', "
                f"chain=<Chain {self.chain.name}>, atoms=[{atoms}], ter={self.ter}")

    def __deepcopy__(self, memo=None):
        if memo is None:
            memo = {}

        if residue := memo.get(id(self)):
            return residue

        residue = replace(self, chain=None)
        memo[id(residue)] = residue
        residue.chain = deepcopy(self.chain, memo)

        return residue

    @property
    def atoms(self) -> list[Atom]:
        if not self.atoms_:
            self.atoms_ = [x for x in self.chain.mol.atoms if x.residue == self]
        return self.atoms_

    @atoms.setter
    def atoms(self, atoms: list[Atom]):
        self.atoms_ = atoms


@dataclass
class Chain:
    idx: int
    name: str
    mol: Molecule | None = None
    residues_: list[Residue] = field(default_factory=list)

    def __repr__(self):
        residues = ', '.join([f'<Residue {x.name}{x.number}>' for x in self.residues])
        return f"Chain(name='{self.name}', residues=[{residues}])"

    def __deepcopy__(self, memo=None):
        if memo is None:
            memo = {}

        if chain := memo.get(id(self)):
            return chain

        chain = replace(self, mol=None)
        memo[id(self)] = chain
        chain.mol = deepcopy(self.mol, memo)

        return chain

    @property
    def residues(self) -> list[Residue]:
        if not self.residues_:
            self.residues_ = [x for x in self.mol.residues if x.chain == self]
        return self.residues_

    @residues.setter
    def residues(self, residues: list[Residue]):
        self.residues_ = residues


class MolType(str, Enum):
    bio = 'BIOPOLYMER'
    protein = 'PROTEIN'
    na = 'NUCLEIC_ACID'
    sugar = 'SACCHARIDE'
    small = 'SMALL'
    unknown = 'UNKNOWN'


@dataclass
class Molecule:
    name: str
    comment: str = ''
    type: MolType = MolType.small
    charge: int = None
    charge_type: str = 'NO_CHARGES'
    data: dict = field(default_factory=dict)
    chains: list[Chain] = field(default_factory=list)
    residues: list[Residue] = field(default_factory=list)
    atoms: list[Atom] = field(default_factory=list)
    bonds: list[Bond] = field(default_factory=list)

    def __add__(self, other: Molecule) -> Molecule:
        mol = self.__deepcopy__()
        other = other.__deepcopy__()

        # Set molecular properties
        mol.name = '-'.join([self.name, other.name])
        mol.comment = ''
        mol.type = MolType.small
        if other.charge is None:
            mol.charge = None
        if other.charge_type == 'NO_CHARGES':
            mol.charge_type = 'NO_CHARGES'
        mol.data |= other.data

        # Make chains
        chain_names = sorted(set(x.name for x in mol.chains + other.chains))
        mol.chains = [Chain(idx=i, name=name, mol=mol) for i, name in
                      enumerate(chain_names)]

        # Make residues
        mol.residues[-1].ter = other.residues[-1].ter = True
        mol.residues.extend(other.residues)

        for i, residue in enumerate(mol.residues):
            residue.idx = i
            chain = next(x for x in mol.chains if x.name == residue.chain.name)
            residue.chain = chain

        # Make atoms
        mol.atoms.extend(other.atoms)
        for i, atom in enumerate(mol.atoms):
            atom.idx = i
            atom.number = i + 1

        # Make bonds
        mol.bonds.extend(other.bonds)
        mol.bonds.sort(key=lambda x: (x.atom.idx, x.partner.idx))

        return mol

    def __deepcopy__(self, memo=None) -> Molecule:
        if memo is None:
            memo = {}

        if mol := memo.get(id(self)):
            return mol

        mol = replace(self, chains=[], residues=[], atoms=[], bonds=[])
        memo[id(self)] = mol

        mol.chains = [deepcopy(x, memo=memo) for x in self.chains]
        mol.residues = [deepcopy(x, memo=memo) for x in self.residues]
        mol.atoms = [deepcopy(x, memo=memo) for x in self.atoms]
        mol.bonds = [deepcopy(x, memo=memo) for x in self.bonds]

        return mol

    def __repr__(self):
        chains = ', '.join([f'<Chain {x.name}>' for x in self.chains])
        residues = ', '.join([f'<Residue {x.name}{x.number}>' for x in self.residues])
        atoms = ', '.join([f'<Atom {x.number}.{x.name}>' for x in self.atoms])
        bonds = ', '.join([x.__repr__() for x in self.bonds])
        return (f"Molecule(name='{self.name}', "
                f"comment='{self.comment}', "
                f"type='{self.type}', "
                f"charge_type='{self.charge_type}', "
                f"data={self.data}, "
                f"chains=[{chains}], "
                f"residues=[{residues}], "
                f"atoms=[{atoms}], "
                f"bonds=[{bonds}]")


def guess_residue_type(name: str) -> ResType:
    for typename, namelist in RESTYPES.items():
        if name in namelist:
            return ResType(typename)
    return ResType.unknown


def guess_mol_type(mol: Molecule) -> MolType:
    res_types = set(x.type for x in mol.residues)
    if ResType.bio in res_types:
        return MolType.bio
    elif ResType.small in res_types:
        if res_types & {ResType.protein, ResType.dna, ResType.rna}:
            return MolType.bio
        return MolType.small
    elif ResType.protein in res_types:
        return MolType.protein
    elif ResType.dna in res_types or ResType.rna in res_types:
        return MolType.na
    elif ResType.sugar in res_types:
        return MolType.sugar
    else:
        return MolType.small
