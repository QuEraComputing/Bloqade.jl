import platform
from typing import List
import juliacall
from numpy import typename
pyconvert_ptr = juliacall.Main.seval('pyconvert')

class JLType(object):

    def __init__(self, jlname : str, typevars : List = []) -> None:
        self.jlname = jlname
        self.typevars = typevars

        eval_str = self.jlname
        if self.typevars:
            typevars_str = ','.join(map(lambda x: x.jlname, self.typevars))
            eval_str += '{' + typevars_str + '}'
        self.jltype_object = juliacall.Main.seval(eval_str)

    @property
    def typename(self) -> str:
        name = 'JL' + self.jlname
        if self.typevars:
            typevars = ','.join(map(lambda x: x.typename, self.typevars))
            self._typename = name + '[' + typevars + ']'
        else:
            self._typename = name
        return self._typename

    def __repr__(self) -> str:
        return self.typename

    def __getitem__(self, *typevars):
        if self.typevars: # not empty
            raise ValueError('expect UnionAll, got JL' + \
                self.jlname + '[' + self.typevars[0].typename + ']')
        return JLType(self.jlname, typevars)

    def convert(self, pyobj):
        return pyconvert_ptr(self.jltype_object, pyobj)


for i in [8, 16, 32, 64, 128]:
    jlname = 'Int' + str(i)
    pyname = 'JL' + jlname
    globals()[pyname] = JLType(jlname)

if platform.machine().endswith('64'):
    JLInt = JLInt64
else:
    JLInt = JLInt32

for i in [16, 32, 64]:
    jlname = 'Float' + str(i)
    pyname = 'JL' + jlname
    globals()[pyname] = JLType(jlname)

JLVector = JLType('Vector')
JLMatrix = JLType('Matrix')
JLArray = JLType('Array')
