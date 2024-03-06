import re
from typing import Iterable

from .core import MolType
from .mol import *

PATTERN_ATOM = re.compile(
    r'([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)\s+(\w+)')
# PATTERN_BOND = re.compile(r'(\d+)\s+(\d+)\s+(\d+)')
PATTERN_DATA = re.compile(r'> {2}<(.*)>.*\n(.*)\n+')
PATTERN_CHG = re.compile(r'M {2}CHG\s*\d+\s+(\d+)\s+(-?\d+)')


class Reader:
    def __init__(self, name: str):
        self.name = name

    @staticmethod
    def _parse_atom(mol: Molecule, line: str, idx: int = 0):
        x, y, z = float(line[:10]), float(line[10:20]), float(line[20:30])
        element = line[30:33].strip()
        mol.atoms.append(Atom(idx=idx,
                              number=idx + 1,
                              name=element,
                              x=float(x),
                              y=float(y),
                              z=float(z),
                              element=element,
                              type=element,
                              is_het=True,
                              residue=mol.residues[0],
                              ))

    @staticmethod
    def _parse_bonds(mol: Molecule, num_bonds: int, block: list[str]):
        bonds = []
        for i in range(num_bonds):
            line = block[i]
            id1, id2, bond_order = int(line[:3]), int(line[3:6]), int(line[6:9])
            if id1 > id2:
                id1, id2 = id2, id1
            bonds.append((id1, id2, BondType.get_bond_type_by_count(bond_order)))
        parse_bonds(mol, sorted(bonds))

    @staticmethod
    def _parse_charges(mol: Molecule, lines: list[str]):
        for line in lines:
            if line.startswith('M  CHG'):
                idx, charge = PATTERN_CHG.findall(line)[0]
                if charge.startswith('-'):
                    charge = -int(charge[1:])
                else:
                    charge = int(charge)
                mol.atoms[int(idx) - 1].charge = charge
                mol.charge = mol.charge or 0
                mol.charge += charge

    def _get_blocks(self, fio: TextIO) -> Iterable[list[str]]:
        lines = []
        for line in fio:
            lines.append(line)
            if line.startswith('$$$$'):
                yield lines
                lines = []

        yield lines

    def _read_mol(self, block: list[str]) -> Molecule | None:
        name = block[0].strip()
        # software = block[1].strip()
        comment = block[2].strip()
        num_atoms, num_bonds = int(block[3][:3]), int(block[3][3:6])

        mol = Molecule(name=name or self.name,
                       comment=comment,
                       type=MolType.small,
                       )
        mol.chains = [Chain(idx=0, name='A', mol=mol)]
        mol.residues = [Residue(idx=0,
                                name='LIG',
                                number=1,
                                type=ResType.small,
                                chain=mol.chains[0])]

        # Atom block
        for i in range(num_atoms):
            self._parse_atom(mol, block[i + 4], i)

        # Bonds
        self._parse_bonds(mol, num_bonds, block[num_atoms + 4:])

        # Properties
        self._parse_charges(mol, block[num_atoms + num_bonds + 4:])

        # Data fields
        data_block = ''.join(block[num_atoms + num_bonds + 4:])
        mol.data = dict(PATTERN_DATA.findall(data_block))

        return mol

    def read(self, fio: TextIO) -> Iterable[Molecule]:
        for block in self._get_blocks(fio):
            if not block:
                continue

            mol = self._read_mol(block)
            if mol:
                yield mol


class Writer:
    def __init__(self):
        self.mol = None

    def _header_block(self) -> list[str]:
        return [self.mol.name,
                ' YinfoCloud',
                self.mol.comment,
                f'{len(self.mol.atoms):>3}{len(self.mol.bonds):>3}  0  0  0  0  0  0  0  0999 V2000']

    def _atoms_block(self) -> list[str]:
        charge_mapping = {
            0: 0,
            3: 1,
            2: 2,
            1: 3,
            # 'doublet radical': 4,
            -1: 5,
            -2: 6,
            -3: 7
        }
        lines = []
        for atom in self.mol.atoms:
            lines.append(
                f'{atom.x:>10.4f}{atom.y:>10.4f}{atom.z:>10.4f} {atom.element:<3} '
                f'0  {charge_mapping.get(int(atom.charge or 0), 0)}  0  0  0  0')
        return lines

    def _bonds_block(self) -> list[str]:
        bond_type_mapping = {
            BondType.single: 1,
            BondType.double: 2,
            BondType.triple: 3,
            BondType.aromatic: 4,
            BondType.unknown: 8
        }
        lines = []
        for bond in self.mol.bonds:
            lines.append(f'{bond.atom.number:>3}'
                         f'{bond.partner.number:>3}'
                         f'{bond_type_mapping.get(bond.type, 1):>3}'
                         '  0  0  0  0')
        return lines

    def _properties_block(self) -> list[str]:
        lines = []
        # i = 0
        # for atom in self.mol.atoms:
        #     if atom.charge is not None and int(atom.charge) != 0:
        #         i += 1
        #         lines.append(f'M  CHG{i:>3}{atom.number:>4}{int(atom.charge):>4}')
        lines.append('M  END')
        return lines

    def _data_block(self) -> list[str]:
        lines = []
        for key, value in self.mol.data.items():
            lines.append(f'>  <{key}>')
            lines.append(value)
            lines.append('')
        lines.append('$$$$')
        return lines

    def write(self, mol: Molecule) -> str:
        self.mol = mol
        lines = self._header_block()
        lines += self._atoms_block()
        lines += self._bonds_block()
        lines += self._properties_block()
        lines += self._data_block()
        return '\n'.join(lines) + '\n'
