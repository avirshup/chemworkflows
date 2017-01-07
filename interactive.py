class UserInteraction(object):
    def __init__(self):
        self.__name__ = self.__class__.__name__


class SelectAtomsFromOptions(UserInteraction):
    def __call__(self, pdbfile, choices):
        pass


class ProteinMinimizationDisplay(UserInteraction):
    def __call__(self, initial_energy,
                 final_energy,
                 mm_minimization,
                 trajectory,
                 rmsd):
        pass
