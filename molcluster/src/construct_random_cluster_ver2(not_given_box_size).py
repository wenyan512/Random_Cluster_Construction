import subprocess
import time
from pathlib import Path
from pydantic import BaseModel
from typing import Literal
from bcmol.mol import read_file
import math
import numpy as np
from multiprocessing import Pool

TEMPLATE = """\
tolerance 2.0
filetype pdb
output model.pdb

"""

TEMPLATE_COMPONENT = """\
structure {file}
    number {count}
    inside {shape} {size}
end structure

"""


# pydantic
class Component(BaseModel):
    file: str
    count: int


class Input(BaseModel):
    shape: Literal["sphere", "box"] = "box"
    components: list[Component]
    size: tuple[int, int, int] | None = None

    @property
    def total_count(self):
        return sum(component.count for component in self.components)


class Output(BaseModel):
    file: str
    size: tuple[int, int, int]


def prepare_components(components: list[Component]):
    for component in components:
        file = Path(component.file)
        mol = next(read_file(file))
        mol.write(file.with_suffix(".pdb"))


def config(components: list[Component], size: tuple[int, int, int] = None, shape: Literal["sphere", "box"] = "box"):
    if shape == "box":
        size = f"0. 0. 0. {size[0]} {size[1]} {size[2]}"
    else:
        shape = "sphere"
        size = f"0. 0. 0. {size[0] / 2:d}"

    content = TEMPLATE
    for component in components:
        content += TEMPLATE_COMPONENT.format(
            file=Path(component.file).with_suffix(".pdb"),
            count=component.count,
            shape=shape,
            size=size,
        )

    Path("config.inp").write_text(content)


def run_packmol() -> bool:
    process = subprocess.Popen("packmol < config.inp", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    try:
        process.communicate(timeout=10)
        return True
    except Exception as e:
        process.kill()
        # print(e)
        return False


def find_initial_sizes(components: list[Component], total_num: int) -> tuple:
    max_mol_size = 0
    min_mol_size = float('inf')

    for component in components:
        mol_info = Path(component.file.replace('.mol2', '.pdb')).read_text()
        mol_xyz = np.array([list(map(float, i.split()[6:9])) for i in mol_info.split('\n') if i.startswith('HETATM')])

        mol_size = np.amax(mol_xyz, axis=0) - np.amin(mol_xyz, axis=0)
        mol_size = np.max(mol_size)

        max_mol_size = max(max_mol_size, mol_size)  # update the cube box size for largest mol
        min_mol_size = min(min_mol_size, mol_size)  # update the cube box size for smallest mol

    # largest mol in one row
    initial_max_size = math.ceil(max_mol_size * total_num)
    # smallest mol in the cube box
    initial_min_size = math.floor(min_mol_size * (total_num ** (1 / 3)))

    return initial_min_size, initial_max_size


def packmol_worker(args: tuple[tuple[int, int, int], list[Component]]) -> bool:
    size, components = args
    config(components, shape="box", size=size)
    return run_packmol()


def find_optimal_size(components: list[Component], initial_sizes: tuple, initial_step: int = 10) -> int:
    min_size, max_size = initial_sizes  # 30, 50
    step = initial_step  # 5

    while step >= 1:
        if step > (max_size - min_size):
            step = 1
        sizes = [i for i in range(min_size, max_size + 1, step)]
        args_list = [([size] * 3, components) for size in sizes]

        pool = Pool(processes=1)
        time.sleep(1)
        results = pool.map(packmol_worker, args_list)
        pool.close()
        pool.join()

        success_inds = [ind for ind, success in enumerate(results) if success]

        # all fail
        if len(success_inds) == 0:
            return max_size
        else:
            if success_inds[0] == 0:
                max_size = min_size
                break

            max_size = sizes[success_inds[0]]
            min_size = sizes[success_inds[0] - 1]

            step //= 2

    return max_size


def main(inp: Input) -> Output:
    prepare_components(inp.components)
    if inp.size == None:
        initial_sizes = find_initial_sizes(components=inp.components, total_num=inp.total_count)
        optimal_size = find_optimal_size(components=inp.components, initial_sizes=initial_sizes, initial_step=100)
        optimal_size = tuple([optimal_size] * 3)
    else:
        optimal_size = inp.size
    config(components=inp.components, shape=inp.shape, size=optimal_size)
    run_packmol()
    return Output(file="model.pdb", size=optimal_size)


if __name__ == "__main__":
    import os

    os.chdir("../test")
    os.environ["PATH"] += ':src/bin/packmol'

    print(main(inp=Input(shape="box",
                         components=[{"file": "methanol.mol2", "count": 1000},
                                     {"file": "water.mol2", "count": 1000}],
                         size=None)))

    # print(main(inp=Input(shape="box",
    #                      components=[{"file": "methanol.mol2", "count": 1000},
    #                                  {"file": "water.mol2", "count": 1000}],
    #                      size=(50, 50, 50))))
