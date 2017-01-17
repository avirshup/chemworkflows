class UserInteraction(object):
    FRONTENDID = None  # str: metadata for creating the MST frontend UI

    def __init__(self):
        self.__name__ = self.__class__.__name__


# Might want to see
# https://gist.github.com/maartenbreddels/3378e8257bf0ee18cfcbdacce6e6a77e
# for async widgets

class SelectAtomsFromOptions(UserInteraction):
    FRONTENDID = None  # str: metadata for creating the MST frontend UI

    def __call__(self, pdbfile, choices):
        import moldesign as mdt
        import ipywidgets as ipy

        # TODO: currently MOCKED
        #m = mdt.read(pdbfile)
        #m.draw(display=True)

        k = choices.keys()[0]

        return {'atom_ids': choices[k]}


class ProteinMinimizationDisplay(UserInteraction):
    def __call__(self, initial_energy,
                 final_energy,
                 trajectory,
                 rmsd):
        print 'Initial:', initial_energy
        print 'Final:', final_energy
        print 'RMSD: ', rmsd
        try:
            trajectory.draw(display=True)
        except:
            pass
