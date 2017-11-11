import numpy as np
import gint


def to_atom_numbers(atoms):
    return [gint.element.find(at).atom_number for at in atoms]


def distribute_electrons(n_electrons, multiplicity):
    """
    Return a tuple of (n_alpha, n_beta) resulting in the multiplicity
    and the electron count provided.
    If this cannot be achieved, raises a ValueError.
    """
    if multiplicity <= 0:
        raise ValueError("The multiplicity needs to be a positive number.")

    if n_electrons < 0:
        raise ValueError("The electron count needs to be zero or larger.")

    spin_twice = multiplicity - 1
    if spin_twice % 2 != n_electrons % 2:
        raise ValueError("Only a system with an even number of electrons can "
                         "have an odd multiplicity and vice versa. This system "
                         "has " + str(n_electrons) + " electrons, but a "
                         "multiplicity of " + str(multiplicity) + " is desired.")

    if spin_twice > n_electrons:
        raise ValueError("A system with " + str(n_electrons) + " electrons cannot "
                         "have a multiplicity larger than " + str(n_electrons + 1) +
                         ". You requested a multiplicity of " + str(multiplicity) +
                         ", however.")

    n_alpha = (n_electrons - spin_twice) // 2 + spin_twice
    n_beta = n_electrons - n_alpha

    assert n_alpha >= 0 and n_beta >= 0
    assert n_alpha - n_beta == spin_twice
    return n_alpha, n_beta


class System:
    """
    Class representing a molecular system, i.e. a structure and a number of electrons.
    """
    def __init__(self, atoms=[], coords=[], electrons=None):
        """
        Initialise a molecular system. By default it is empty (i.e. contains no atoms
        or electrons). In order to construct a non-empty system, specify one or more
        of the optional parameters. Certain combinations are not valid and will
        raise a ValueError or TypeError.

        atoms:   List of all atoms of the structure. Either symbol or atomic numbers
                 or names are accepted.
        coords:  List of the coordinates of the involved atoms.

        electrons       The number of electrons in the system.
                        Can be given as a tuple (n_alpha, n_beta) or as a single
                        value. Then we take n_alpha = electrons // 2 and
                        n_beta = electrons - n_alpha
        multiplicity    The multiplicity (2S+1) of the ground state to compute.
                        By default 1 is chosen for even electron systems and 2
                        for odd electron systems.
        charge          Computed from the electrons and the nuclei in the system
                        or set to 0 (neutral atom).
        """
        if isinstance(atoms, (int, str)):
            atoms = [atoms]
        elif not isinstance(atoms, (list, tuple, np.ndarray)):
            raise TypeError("atoms needs to be a list or tuple")

        if len(atoms) > 0 and len(coords) == 0:
            if len(atoms) == 1:
                coords = [[0, 0, 0]]
            else:
                raise ValueError("coords needs to be specified if more than one atom is "
                                 "in the structure.")

        if len(atoms) != len(coords):
            raise ValueError("Number of atoms and number of coordinates does not agree.")

        self.atom_numbers = np.array(to_atom_numbers(atoms))
        self.coords = np.array(coords)

        if len(coords) > 0 and self.coords.shape[1] != 3:
            raise ValueError("The coords list needs to have exactly 3 items per "
                             "coordinate.")

        if isinstance(electrons, (tuple, list)):
            self.n_alpha, self.n_beta = electrons
            if self.n_alpha < self.n_beta:
                raise ValueError("molsturm assumes that n_alpha >= n_beta at many "
                                 "places. You instructed to construct a system with "
                                 "n_beta (== " + str(self.n_beta) + ") > n_alpha ( == " +
                                 str(self.n_alpha) + "), which is not supported.")
            return

        if electrons is None:
            # Set electron count such that charge cancels to zero
            nuclear_charge = self.atom_numbers.sum()
            electrons = nuclear_charge
        elif not isinstance(electrons, int):
            raise TypeError("electrons needs to be a tuple, list or int")

        # Distribute to alpha and beta with alpha getting more electrons:
        self.n_beta = electrons // 2
        self.n_alpha = electrons - self.n_beta

        assert (electrons % 2 == 0 and self.multiplicity == 1) or \
            (electrons % 2 == 1 and self.multiplicity == 2)

    def adjust_electrons(self, charge=None, multiplicity=None,
                         allow_multiplity_change=True):
        """
        Distribute a certain number of electrons into alpha and beta electrons,
        such that a particular charge and multiplicity is achieved.

        charge         The total charge to achieve in the molecular system after the
                       distribution of electrons. If this parameter is absent,
                       the present value (by default zero) is not altered.
        multiplicity   The multiplicity to achieve after the distribution.
                       If this value is not set, the function will try to retain
                       the original multiplicity. If this cannot be done,
                       then the multiplicity will be 1 for even electron systems
                       and 2 for odd electron systems after the call.
        allow_multiplity_change   If this is set to false, than changing the charge
                       in a way that the current multiplicity value becomes
                       invalid is not allowed and yields a ValueError.
        """
        if charge is None:
            n_elec_count = self.n_electrons
        else:
            nuclear_charge = self.atom_numbers.sum()
            n_elec_count = nuclear_charge - charge
            if n_elec_count < 0:
                raise ValueError("Charge cannot be more positive than the "
                                 "total nuclear charge")

        if multiplicity is None:
            # Try first with original multiplicity.
            try:
                self.n_alpha, self.n_beta = distribute_electrons(n_elec_count,
                                                                 self.multiplicity)
                return
            except ValueError:
                if not allow_multiplity_change:
                    raise ValueError("Changing the charge in the way requested would "
                                     "require a change in multiplicity as well.")

            # Set to 1 for even-electron count and 2 for odd:
            multiplicity = 1 if n_elec_count % 2 == 0 else 2

        self.n_alpha, self.n_beta = distribute_electrons(n_elec_count, multiplicity)

    @property
    def n_atoms(self):
        return len(self.atom_numbers)

    @property
    def multiplicity(self):
        """Return the multiplicity of the system"""
        return self.n_alpha - self.n_beta + 1

    @property
    def is_closed_shell(self):
        """Is the system closed shell"""
        return self.n_alpha == self.n_beta

    @property
    def n_electrons(self):
        return self.n_alpha + self.n_beta

    @property
    def total_charge(self):
        """Return the total resulting charge of the system"""
        return self.atom_numbers.sum() - self.n_electrons

    @property
    def charge(self):
        """Return the total resulting charge of the system"""
        return self.total_charge
