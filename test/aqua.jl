using PorousMaterials
import Aqua

# ambiguity testing finds many "problems" outside the scope of this package
ambiguities=false

# to skip when checking for stale dependencies and missing compat entries
#   Aqua is added in a separate CI job, so (ironically) does not work w/ itself
stale_deps = (ignore=[:Aqua],)

Aqua.test_all(
    PorousMaterials;
    ambiguities=ambiguities, 
    stale_deps=stale_deps
)