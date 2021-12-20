# file with example command line parameters
# more info, see the function init_params() from the *cpp files
# parameters:
#   pLeave: probability predator leaves
#   pArrive: probability predator arrives
#   pAttack: probability predator attacks if it is there
#   alpha: how mortality function decreases with increasing hormone levels (concave, convex, linear)
#   Kmort: how strongly mortality increases with increasing damage
#   Kfec: how strongly fecundity decreases with increasing damage

# executable name       pLeave  pArrive pAttack alpha   Kmort   Kfec
./stress_damage_lh.exe  0.95    0.05    0.5     1.0     0.001   0.01 
