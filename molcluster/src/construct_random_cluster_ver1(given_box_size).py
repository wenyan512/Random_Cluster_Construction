import subprocess
from pathlib import Path
from pydantic import BaseModel
from typing import Literal
from bcmol.mol import read_file

TEMPLATE = '''\
tolerance 2.0
filetype pdb
output model.pdb

'''

TEMPLATE_COMPONENT = '''\
structure {file}
    number {count}
    inside {shape} {size}
end structure

'''


class Component(BaseModel):
    file: str
    count: int


class Input(BaseModel):
    shape: Literal['sphere', 'box'] = 'box'
    components: list[Component]
    size: tuple[int, int, int] | None = None


class Output(BaseModel):
    file: str
    size: tuple[int, int, int]


def prepare_components(components: list[Component]):
    for component in components:
        file = Path(component.file)
        mol = next(read_file(file))
        mol.write(file.with_suffix('.pdb'))


def config(components: list[Component], size: tuple[int, int, int] = None, shape: Literal['sphere', 'box'] = 'box'):
    if shape == 'box':
        size = f'0. 0. 0. {size[0]} {size[1]} {size[2]}'
    else:
        shape = 'sphere'
        size = f'0. 0. 0. {size[0] / 2:d}'

    content = TEMPLATE

    for component in components:
        content += TEMPLATE_COMPONENT.format(file=Path(component.file).with_suffix('.pdb'),
                                             count=component.count,
                                             shape=shape,
                                             size=size)

    Path('./config.inp').write_text(content)


def run_packmol():
    s, o = subprocess.getstatusoutput('packmol < config.inp')
    if s:
        raise RuntimeError(o)


def main(inp: Input) -> Output:
    prepare_components(inp.components)
    config(inp.components, inp.size, inp.shape)
    run_packmol()
    return Output(file='model.pdb', size=inp.size)


if __name__ == '__main__':
    import os

    os.chdir('../test')
    os.environ['PATH'] += ':/home/zhangzy/PycharmProjects/yinfo-tools/molcluster/src/bin'
    print(main(inp=Input(shape='box',
                         components=[{'file': 'methanol.mol2', 'count': 1000},
                                     {'file': 'water.mol2', 'count': 1000}],
                         size=(50, 50, 50))))
