pyjl_block_repr(self) = repr(self)

function init_jlwrap_block()
    pybuiltins.exec(pybuiltins.compile("""
    class YaoBlock(juliacall.AnyValue):

    def __repr__(self) -> str:
        return self._jl_callmethod($(pyjl_methodnum(pyjl_block_repr)))
    """)
end
