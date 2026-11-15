# Per-node reduced-cost arc fixing.
#
# removedArcs (global BitVector, read by populateReducedCosts: a removed arc gets Inf
# reduced cost, so the C++ labeling never uses it) holds the arcs removed at the
# CURRENT node. A Coluna record unit stores a copy with every node and restores it
# before the node (or a strong-branching child) is processed, as MasterCutsUnit does
# for the cuts (pattern: Coluna test improve_relax_callback.jl, ToyNodeInfoUnit).
# Arcs fixed at the root are in the root's record, so every node inherits them.
#
# The fixing runs in CutCleanup (BeforeCutGenAlgo: after every column generation, at
# the root and in every strong-branching phase-2 child), with that node's LP value.
# Before (until 2026-11-15) it ran only at the root, in the cut callback.

const ClB = Coluna.ColunaBase
const ClA = Coluna.Algorithm

struct ArcFixUnit <: ClB.AbstractRecordUnit end
ClB.storage_unit(::Type{ArcFixUnit}, _) = ArcFixUnit()

# Outside the root, fix only when the node's LP value rose by at least NodeFixMinRise
# since the last fixing in this branch: each fixing is a full exact labeling (~0.09 s on
# E-n101) and between the cut rounds of a strong-branching child the LP value barely
# moves (fixing at every call: 324 calls, 29.6 s in a 10-min E-n101 run). The root fixes
# at every call (few calls, and its arcs are removed in the whole tree).
const NodeFixMinRise = 0.5
const lastFixZLP     = Ref(-Inf)            # LP value of the last fixing (current branch)

struct ArcFixRecord <: ClB.AbstractRecord
    removed::BitVector                      # copy, never the live vector
    lastFixZLP::Float64
end

struct ArcFixKey <: ClA.AbstractStorageUnitKey end
ClA.key_from_storage_unit_type(::Type{ArcFixUnit}) = ArcFixKey()
ClA.record_type_from_key(::ArcFixKey) = ArcFixRecord
ClB.record_type(::Type{ArcFixUnit}) = ArcFixRecord
ClB.storage_unit_type(::Type{ArcFixRecord}) = ArcFixUnit

ClB.record(::Type{ArcFixRecord}, id::Int, form::Coluna.MathProg.Formulation, ::ArcFixUnit) =
    ArcFixRecord(copy(removedArcs), lastFixZLP[])

function ClB.restore_from_record!(form::Coluna.MathProg.Formulation, ::ArcFixUnit,
                                  record::ArcFixRecord)
    copyto!(removedArcs, record.removed)
    lastFixZLP[] = record.lastFixZLP
    return
end

# The unit is declared by CutCleanup, which runs inside BeforeCutGenAlgo. Coluna only
# restores the units of algorithms that are not "managers" (collect_units_to_restore!,
# Algorithm/interface.jl), and every conquer algorithm is a manager
# (ismanager(::AbstractConquerAlgorithm) = true), BeforeCutGenAlgo included: without
# this line the unit is recorded with every node but never restored (each node then
# started from the arcs fixed at the previous one: invalid fixing, E-n51 lost its
# optimum). Coluna only checks ismanager of the wrapped algorithm (CutCleanup) when
# running it (conquer.jl, run_before_cutgen_user_algo!), so this changes nothing else.
Coluna.Algorithm.ismanager(::Coluna.Algorithm.BeforeCutGenAlgo) = false

# Statistics (count(removedArcs) at the end is only the last node's set).
const nNodeFixCalls   = Ref(0)
const nNodeFixRemoved = Ref(0)   # arcs removed, summed over the calls
const maxNodeRemoved  = Ref(0)   # largest set of removed arcs at a node
const nNodeFixSkipped = Ref(0)   # calls skipped (LP value did not rise enough)

# Called from CutCleanup.run! with the node's conquer output.
function nodeArcFixing!(conquer_output)
    ArcFixing || return
    UB  = ubArcFixing[]
    zLP = Coluna.ColunaBase.getvalue(ClA.get_lp_primal_bound(conquer_output))
    # Not finite, or the node will be pruned (zLP >= UB): a negative gap would remove
    # every arc, so do nothing.
    (isfinite(UB) && isfinite(zLP) && zLP < UB) || return
    if rootDone[] && zLP < lastFixZLP[] + NodeFixMinRise
        nNodeFixSkipped[] += 1
        return
    end
    newly = fixArcs(zLP, inst)
    lastFixZLP[] = zLP
    nNodeFixCalls[]   += 1
    nNodeFixRemoved[] += newly
    maxNodeRemoved[]   = max(maxNodeRemoved[], count(removedArcs))
    return
end

nodeArcFixingReset!() = (nNodeFixCalls[] = 0; nNodeFixRemoved[] = 0; maxNodeRemoved[] = 0;
                         nNodeFixSkipped[] = 0; lastFixZLP[] = -Inf; nothing)
nodeArcFixingPrint() = println("Node arc fixing: ", nNodeFixCalls[], " calls (", nNodeFixSkipped[],
    " skipped), ",
    nNodeFixRemoved[], " arcs removed in total, at most ", maxNodeRemoved[], " at a node")
