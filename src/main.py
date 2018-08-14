#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# =============================================================================
# Begins main.py
#
# Main program that performs 3D Monte Carlo simulations on lattice proteins.
# =============================================================================

import math, random
from src.energy import energy
from src.changeit import changeit

def main(T, amp, period, iteration, protein_name, random_conf = False):
    """
    Main program that performs 3D Monte Carlo simulations of lattice proteins.

    The call of main function is:
        main(T, amp, period, iteration, protein_name),
    and one optional flag "random_conf". `True` generates a random \
    self-avolding initial conformation, and `False` uses all straight \
    as initial conformation. By default, `random_conf = False`.

    It returns a tuble of five lists:
        (acceptance, etotal, conf_list, native, all_contact).
    """

# =============================================================================
#   1.  Input check.
# =============================================================================

    assert isinstance(iteration, int), "Iteration Number is not an integer!"
    assert T > 0, "Temperature must be greater than 0!"

# =============================================================================
#   2. Initialization.
# =============================================================================

    if protein_name == "mer4":
        from proteins.mer4 import seq as seq
        from proteins.mer4 import native_list as native_list
    elif protein_name == "mer48A":
        from proteins.mer48A import seq as seq
        from proteins.mer48A import native_list as native_list
    elif protein_name == "mer27":
        from proteins.mer27 import seq as seq
        from proteins.mer27 import native_list as native_list
    elif protein_name == "mer15":
        from proteins.mer15 import seq as seq
        from proteins.mer15 import native_list as native_list
    elif protein_name == "chignolin":
        from proteins.chignolin import seq as seq
        from proteins.chignolin import native_list as native_list
    else: raise AssertionError("Invalid protein name. Please check the input!")

    # Generate initial conformation.
    if random_conf == False:
        conformation = "A" * (len(seq) - 1)
    elif random_conf == True:
        # Placeholder.
        conformation = "A" * (len(seq) - 1)
    else: raise AssertionError("Invalid random_conf argument!")
    # Length check
    assert len(seq) == len(conformation) + 1

    # A Boolean list stores decisions.
    acceptance = []
    # A float list stores energies.
    etotal = []
    # Two integer lists stores number of contacts.
    num_native_list = []
    num_all_list = []
    # A string list stores conformations.
    conf_list = []

    # Temporary local variables.
    e = 0
    contact = 0
    native_contact = 0
    metro = 0

    # A percentage marker.
    if iteration > 100:
        percentage = iteration//100
    else:
        percentage = 1

# =============================================================================
#   3. Perform simulation.
# =============================================================================

    for i in range(iteration):

        # A percentage output to screen
        if i % percentage == 0:
            print (i/percentage, "% done.")
        if (i//period)%2 == 0:
            temp = T + amp
        else:
            temp = T - amp


        # 3.1 Propose a new conformation
        newconformation = changeit(conformation)

        # 3.2 Compute the contact energy of new conformation.
        (newe, newcontact, newnative_contact) = energy(newconformation, seq,
                                                       native_list)

        # 3.3 Proposed move decisions based on Metropolis probability.

        # Corrected on 06/28/18 commit.
        # Self-avoiding check was moved into `changeit` module.
        # Reject the move if the new conformation is not self-avoiding.
        # if newe == float('nan'):
        #     acceptance.append(0)
        # Accept the move if energy drops or unchanged.
        if newe <= e:
            conformation = newconformation
            e = newe
            contact = newcontact
            native_contact = newnative_contact
#            print("Accept")
            acceptance.append(1)
        # Accept the move at the Metropolis probability if energy increases.
        # According to Miyazawa-Jernigen (1996), the room temperature (300 K) was
        #   assumed when converting the RT energy units to kcal/mol.
        else:
            beta = 300.0/temp
            metro = math.exp(beta * (e - newe))
            if random.random() < metro:
                conformation = newconformation
                e = newe
                contact = newcontact
                native_contact = newnative_contact
#                print("Accept")
                acceptance.append(1)
            # Reject the move otherwise.
            else: acceptance.append(0)

        etotal.append(e)
        num_native_list.append(native_contact)
        num_all_list.append(contact)
        conf_list.append(conformation)

        # print (i, newconformation, metro)

# =============================================================================
#   4. Return conformational list and other trajectory information.
# =============================================================================
    return acceptance, etotal, conf_list, num_native_list, num_all_list