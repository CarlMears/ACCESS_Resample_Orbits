class Computational_Library:
    def __init__(self, dll_path: str) -> None:
        self.dll = None
        
        try:
            find_problem_dll(dll_path)
            self.dll = WinDLL(dll_path)
        except OSError:
            raise RuntimeError(
                f"Did not find computational library DLL. Searched: {dll_path}"
            )
        except:
            raise

        # self.set_hook = self.dll.set_table
        # self.set_hook.argtypes = [POINTER(c_int)]
        # self.set_hook.restype = c_int

        # self.add_hook = self.dll.add_table
        # self.add_hook.argtypes = [POINTER(Input), POINTER(Output)]
        # self.add_hook.restype = c_int

        # self.del_hook = self.dll.del_table
        # self.del_hook.argtypes = None
        # self.del_hook.restype = c_int
