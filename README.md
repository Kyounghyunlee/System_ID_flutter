# Flutter rig system identification

Linear system ID code computes Hopf point and Modal vectors of the grey-box model (using unsteady flutter model)
from the free-decay response of the aero foil.
Additionally it computes input form of Julia code that computes Center manifold and Hopf normal form

Nonlinear system ID code computes nonlinear stiffness parameters ka2 and ka3 by minimizing the prediction error
using normal form theory. This result can be used as a initial searching point of function optimize_NLstiff(x1,x2), optimize_NLstiff2(x1,x2), Julia package HopfNormalForm.
