import re
from collections import Counter
from typing import Iterable

from .mol import *

PATTERN_PDB_CHARGE = re.compile(r'([0-9])([+-]?)')


class Reader:
    def __init__(self, name: str = None):
        self.name = name
        self.atoms_dict = {}
        self.residues_dict = {}
        self.chains_dict = {}

    def _parse_atom(self, mol: Molecule, line: str, idx: int = 0):
        """
        References:
            - https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
        """
        # Construct chain
        chain_name = line[21:22].strip()
        chain = self.chains_dict.get(chain_name)
        if chain is None:
            chain = Chain(idx=0, name=chain_name, mol=mol)
            mol.chains.append(chain)
            self.chains_dict[chain_name] = chain

        # Construct residue
        residue_name = line[17:20].strip()
        residue_number = line[22:26]
        key = (residue_name, residue_number, chain.name)
        residue = self.residues_dict.get(key)
        if residue is None:
            residue = Residue(idx=0,
                              name=residue_name,
                              number=int(residue_number),
                              chain=chain,
                              type=guess_residue_type(residue_name))
            mol.residues.append(residue)
            self.residues_dict[key] = residue

        # Construct atom
        charge = line[78:80].strip()
        if charge:
            charge, sign = re.findall(PATTERN_PDB_CHARGE, charge)[0]
            charge = float(f'{sign}{charge}')
        else:
            charge = None

        atom = Atom(idx=idx,
                    number=int(line[6:11]),
                    name=line[12:16].strip(),
                    alt=line[16:17].strip(),
                    insertion=line[26:27].strip(),
                    x=float(line[30:38]),
                    y=float(line[38:46]),
                    z=float(line[46:54]),
                    occupancy=float(line[54:60]),
                    bfactor=float(line[60:66]),
                    segid=line[72:76].strip(),
                    element=line[76:78].strip(),
                    charge=charge,
                    is_het=line[0:6].strip() == 'HETATM',
                    residue=residue,
                    )
        mol.atoms.append(atom)

    @staticmethod
    def _parse_bonds(line: str) -> list[tuple[int, int, BondType]]:
        """
        References:
            - https://www.wwpdb.org/documentation/file-format-content/format33/sect10.html#CONECT
        """
        line = f'{line:31}'
        atom = int(line[6:11].strip())
        partners = [line[11:16].strip(),
                    line[16:21].strip(),
                    line[21:26].strip(),
                    line[26:31].strip()]
        partners = [int(x) for x in filter(None, partners)]
        partners = Counter(partners)

        return [(atom, partner, BondType.get_bond_type_by_count(count))
                for partner, count in partners.items()]

    @staticmethod
    def _get_blocks(fio) -> Iterable[list[str]]:
        lines = []
        for line in fio:
            lines.append(line)
            if line.startswith('MODEL'):
                lines = []
            elif line.startswith('ENDMDL'):
                yield lines
        yield lines

    def _read_mol(self, lines: list[str]) -> Molecule | None:
        self.atoms_dict = {}
        self.residues_dict = {}
        self.chains_dict = {}
        mol = Molecule(name=self.name)

        # Handle atoms
        i = 0
        bonds = []
        for line in lines:
            if line.startswith(('ATOM', 'HETATM')):
                self._parse_atom(mol, line, i)
                i += 1
            elif line.startswith('TER'):
                mol.atoms[-1].residue.ter = True
            elif line.startswith('CONECT'):
                bonds.extend(self._parse_bonds(line))

        if not mol.atoms:
            return

        # Handle chains & residues
        for i, chain in enumerate(mol.chains):
            chain.idx = i
        for i, residue in enumerate(mol.residues):
            residue.idx = i

        # Handle mol type
        mol.type = guess_mol_type(mol)

        # Handle bonds
        if bonds:
            parse_bonds(mol, sorted(bonds))

        return mol

    def read(self, fio: TextIO) -> Iterable[Molecule]:
        for block in self._get_blocks(fio):
            if not block:
                continue

            mol = self._read_mol(block)
            if mol:
                yield mol


class Writer:
    def __init__(self, increase_tercount: bool = True):
        self.mol = None
        self.increase_tercount = increase_tercount

    @staticmethod
    def normal_record(atom: Atom):
        """
        References
            https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
            https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#HETATM
        """
        residue_ = atom.residue
        record_type = 'HETATM' if atom.is_het else 'ATOM'
        atom_name = f'{atom.name:<3}'
        # if (charge := atom_.charge) > 0:
        #     charge = f'{charge}+'
        # elif charge < 0:
        #     charge = f'{abs(charge)}-'
        # else:
        #     charge = ''
        charge = ''
        return (f'{record_type:<6}{atom.number:>5} {atom_name:>4}{atom.alt:1}'
                f'{residue_.name:>3} {residue_.chain.name:1}'
                f'{residue_.number:>4}{atom.insertion:1}   '
                f'{atom.x:>8.3f}{atom.y:>8.3f}{atom.z:>8.3f}'
                f'{atom.occupancy:>6.2f}{atom.bfactor:>6.2f}      '
                f'{atom.segid:<4}{atom.element:>2} {charge:<2}\n')

    def ter_record(self, atom: Atom):
        """
        References
            https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#TER
        """
        residue_ = atom.residue
        atom_number = atom.number + 1 if self.increase_tercount else atom.number
        return (f'TER   {atom_number:>5}      {residue_.name:>3} '
                f'{residue_.chain.name:1}{residue_.number:>4}{atom.insertion:1}\n')

    @staticmethod
    def connect_records(bonds: list[Bond]) -> str:
        content = ''
        for bond in bonds:
            bonding_atoms = ''.join([f'{bond.partner.number:>5}'] * bond.order)
            content += f'CONECT{bond.atom.number:>5}{bonding_atoms}\n'
        return content

    def write(self, mol: Molecule, idx: int = None, write_conects: bool = True) -> str:
        self.mol = mol.clone()

        content, offset = '', 0
        last_atoms = [x.atoms[-1].idx for x in self.mol.residues if x.ter]
        for i, atom in enumerate(self.mol.atoms, 1):
            atom.number = offset + i
            content += self.normal_record(atom)

            # Write coordinate records
            if atom.idx in last_atoms:
                if self.increase_tercount:
                    offset += 1
                content += self.ter_record(atom)

        # Whether to write `MODEL...ENDMDL` records
        if idx is not None:
            content = f'MODEL  {idx:>8}\n' + content + 'ENDMDL\n'

        # Write connectivity records
        if write_conects:
            content += self.connect_records(self.mol.bonds)

        return content
