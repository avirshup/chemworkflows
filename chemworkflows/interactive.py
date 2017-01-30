class UserInteraction(object):
    def __init__(self):
        self.__name__ = self.__class__.__name__


# Might want to see
# https://gist.github.com/maartenbreddels/3378e8257bf0ee18cfcbdacce6e6a77e
# for async widgets


class SelectAtomsFromOptions(UserInteraction):
    def __call__(self, choices):
        if len(choices) > 1:

            print 'Choose your ligand!'
            choicelist = choices.items()

            for i, (choice, value) in enumerate(choicelist):
                print '%d) %s: %s' % (i+1, choice, value)

            userinput = ''
            while not (userinput.isdigit() and 1 <= int(userinput) <= len(choices)):
                userinput = raw_input('Which ligand (1-%d)? ' % len(choices))

            index = int(userinput) - 1
            name, atom_ids = choicelist[index]

        else:
            name = choices.keys()[0]
            atom_ids = choices.values()[0]


        return {'ligandname': name,
                'atom_ids': atom_ids}


