import re

from .mol import *

# atom_id, atom_name, x, y, z, atom_type, residue_idx, residue_name, residue_number, charge
PATTERN_TRIPOS_ATOM_RECORD = re.compile(
    r'(\d+)\s+(\S+)\s+([\d.-]+)\s+([\d.-]+)\s+([\d.-]+)\s+([\w.]+)\s+(\d+)\s+(\S{1,3}|\D+)(\d+)?\s+(-?\d+\.?\d*)')
    # r'(\d+)\s+(\S+)\s+([\d.-]+)\s+([\d.-]+)\s+([\d.-]+)\s+([\w.]+)\s+(\d+)\s+(\S{1,3}|[^0-9]+)(\d+)?\s+(-?\d+\.\d+)')
PATTERN_TRIPOS_BOND_RECORD = re.compile(r'\d+\s+(\d+)\s+(\d+)\s+(\w+)')


class Reader:
    def __init__(self, name: str):
        self.name = name
        self.atoms_dict = {}

    @staticmethod
    def _parse_atom(mol: Molecule, atom: tuple, idx: int = 0):
        """
        Args:
            mol (RawMolecule):
            atom (tuple): atom_id, atom_name, x, y, z, atom_type, residue_number, residue_name, charge
        """
        # Construct residue
        number = int(atom[6])
        residue = Residue(idx=number - 1,
                          name=atom[7],
                          number=int(atom[8] or number),
                          chain=mol.chains[0],
                          )
        residue.type = guess_residue_type(residue.name)
        if residue in mol.residues:
            residue = mol.residues[mol.residues.index(residue)]
        else:
            mol.residues.append(residue)

        atom = Atom(idx=idx,
                    number=int(atom[0]),
                    name=atom[1],
                    x=float(atom[2]),
                    y=float(atom[3]),
                    z=float(atom[4]),
                    element=atom[5].split('.')[0],
                    type=atom[5],
                    charge=float(atom[9]),
                    is_het=True,
                    residue=residue,
                    )
        mol.atoms.append(atom)

    @staticmethod
    def _parse_bonds(mol: Molecule, num_bonds: int, block: list[str]):
        bonds = []
        for i in range(num_bonds):
            id1, id2, bond_type = PATTERN_TRIPOS_BOND_RECORD.findall(block[i])[0]
            id1, id2, bond_type = int(id1), int(id2), BondType(bond_type)
            if id1 > id2:
                id1, id2 = id2, id1
            bonds.append((id1, id2, bond_type))
        parse_bonds(mol, sorted(bonds))

    @staticmethod
    def _get_blocks(fio) -> Iterator[list[str]]:
        seen, block = False, []
        for line in fio:
            if line.startswith('##########'):
                if seen:
                    yield block
                    block = []
                    seen = False
            elif line.startswith('@<TRIPOS>MOLECULE'):
                if seen:
                    yield block
                    block = []
                seen = True
            block.append(line)

        if block:
            yield block

    def _read_mol(self, block: list[str]) -> Molecule:
        """
        References
            https://chemicbook.com/2021/02/20/mol2-file-format-explained-for-beginners-part-2.html

        TODO:
            Use a more efficient way to read mol instead of slicing the block many times.
        """
        mol = Molecule(name=self.name)

        idx = 0
        # UCSF Dock scores saved in comments
        for i, line in enumerate(block):
            if line.startswith('##########'):
                key, value = line[10:].split(':')[:2]
                mol.data[key.strip()] = value.strip()
            elif line.startswith('@<TRIPOS>MOLECULE'):
                idx = i
                break

        # @<TRIPOS>MOLECULE
        block = block[idx:]
        mol.name = block[1].strip() or mol.name
        counts = block[2].split()
        num_atoms, num_bonds = int(counts[0]), int(counts[1])
        mol.type = block[3].strip()
        mol.charge_type = block[4].strip()
        mol.comment = '' if block[5].startswith('@<TRIPOS>ATOM') else block[5].strip()
        mol.chains = [Chain(idx=0, name='A', mol=mol)]

        # @<TRIPOS>ATOM
        idx = next(i for i, line in enumerate(block) if line.startswith('@<TRIPOS>ATOM'))
        block = block[idx + 1:]
        for i in range(num_atoms):
            atom = PATTERN_TRIPOS_ATOM_RECORD.findall(block[i])[0]
            self._parse_atom(mol, atom, i)

        # Handle residues
        for i, residue in enumerate(mol.residues):
            residue.idx = i
        mol.residues[-1].ter = True

        # @<TRIPOS>BOND
        if num_bonds:
            idx = next(i for i, x in enumerate(block) if x.startswith('@<TRIPOS>BOND'))
            block = block[idx + 1:]
            self._parse_bonds(mol, num_bonds, block)

        return mol

    def read(self, fio: TextIO) -> Iterator[Molecule]:
        for block in self._get_blocks(fio):
            if not block:
                continue
            elif block[0] == '\n':
                block = block[1:]

            mol = self._read_mol(block)
            if mol:
                yield mol


class Writer:
    def __init__(self):
        self.mol = None

    def _data_section(self) -> list[str]:
        lines = ['\n']
        for k, v in self.mol.data.items():
            lines.append(f'##########	{k}:\t{v}\n')
        lines.append('\n')
        return lines

    def _mol_section(self) -> list[str]:
        return ['@<TRIPOS>MOLECULE\n',
                f'{self.mol.name}\n',
                f' {len(self.mol.atoms)} {len(self.mol.bonds)} {len(self.mol.residues)} 0 0\n',
                f'{self.mol.type}\n',
                f'{self.mol.charge_type or "NO_CHARGES"}\n',
                f'{self.mol.comment}\n']

    def _atom_section(self) -> list[str]:
        lines = ['@<TRIPOS>ATOM\n']
        for atom in self.mol.atoms:
            atom_name = f'{atom.name:<3}'
            lines.append(f'{atom.number:>7} {atom_name:>4} {" " * 3} '
                         f'{atom.x:>9.4f} {atom.y:>9.4f} {atom.z:>9.4f} '
                         f'{atom.type:<5} {atom.residue.idx + 1:>3} '
                         f'{atom.residue.name:>4}{atom.residue.number:<5} '
                         f'{atom.charge or 0:>9.4f}\n')
        return lines

    def _bond_section(self) -> list[str]:
        lines = ['@<TRIPOS>BOND\n']
        for i, bond in enumerate(self.mol.bonds, 1):
            lines.append(f'{i:>6} {bond.atom.number:>5} {bond.partner.number:>5} {bond.type:>4}\n')
        return lines

    def write(self, mol: Molecule) -> str:
        self.mol = mol
        lines = []
        # Data
        if self.mol.data:
            lines += self._data_section()

        # RawMolecule info (@<TRIPOS>MOLECULE)
        lines += self._mol_section()

        # Coordinates (@<TRIPOS>ATOM)
        lines += self._atom_section()

        # Bonds (@<TRIPOS>BOND)
        if self.mol.bonds:
            lines += self._bond_section()

        return ''.join(lines)
