class YaoBlock(juliacall.AnyValue):

    def __repr__(self) -> str:
        return self._jl_callmethod()
pyjl_methodnum(pyjlany_repr)