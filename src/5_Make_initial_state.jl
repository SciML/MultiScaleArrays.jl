include("./5bis_FunctAndStruct_Make_initial_state.jl")

ahostcell  =  CellHostGenot([1.0 ,AllelesParam.TraitsValues.MaxE_sp1[AllelesParam.AlleleINI.MaxE_sp1[1]]*0.95, 0.0 ]
                            ,:aCellHostGenot
                            ,Int(1)
                            ,typeof(SimParam.Nsp)(1)
                            ,allelesHostSymb[1]
                            ,MakePhenotypesFromAlleles(allelesHostSymb[1],   typeof(SimParam.Nsp)(1) )
                            ,          GetMuEachTrait(allelesHostSymb[1], convert(typeof(SimParam.Nsp),1) , Val(true ) , AllelesParam.MuMu, TraitsValues, AllelesParam.TraitAffectHostsSymbiontsCells)  # MuEachTraitCell
                            ,          GetMuEachTrait(allelesHostSymb[1], convert(typeof(SimParam.Nsp),1) , Val(false) , AllelesParam.MuMu, TraitsValues, AllelesParam.TraitAffectHostsSymbiontsCells)  # MuEachTraitHost
                            ,GetMuTot( GetMuEachTrait(allelesHostSymb[1], convert(typeof(SimParam.Nsp),1) , Val(true ) , AllelesParam.MuMu, TraitsValues, AllelesParam.TraitAffectHostsSymbiontsCells)) # MuTotCell
                            ,GetMuTot( GetMuEachTrait(allelesHostSymb[1], convert(typeof(SimParam.Nsp),1) , Val(false) , AllelesParam.MuMu, TraitsValues, AllelesParam.TraitAffectHostsSymbiontsCells)) # MuTotHostOnly
#                             ,MakeMapMutatingTraits(allelesHostSymb[1], true)
#                             ,MakeMapMutatingTraits(allelesHostSymb[1], false)
                            )

asymbiontS = [
              CellSymbGenot([SimParam.InitialState.NsymbCells ,  AllelesParam.TraitsValues[Symbol(:MaxE_sp,sp)][AllelesParam.AlleleINI[Symbol(:MaxE_sp,sp)]]*0.95 ]
                            ,:aSymbGeno
                            ,Int(1)
                            ,typeof(SimParam.Nsp)(sp)
                            ,allelesHostSymb[sp]
                            ,MakePhenotypesFromAlleles(allelesHostSymb[sp],sp)
                            ,         GetMuEachTrait(allelesHostSymb[sp], sp , Val(true) , AllelesParam.MuMu, TraitsValues, AllelesParam.TraitAffectHostsSymbiontsCells)# MuEachTraitCell # MuEachTrait
                            ,GetMuTot(GetMuEachTrait(allelesHostSymb[sp], sp , Val(true) , AllelesParam.MuMu, TraitsValues, AllelesParam.TraitAffectHostsSymbiontsCells))# MuEachTraitCell # MuEachTrait
#                             ,MakeMapMutatingTraits(allelesHostSymb[sp], true)
                            )
for sp in typeof(SimParam.Nsp)(2):typeof(SimParam.Nsp)(SimParam.Nsp)]

acellgenoS = [ahostcell,asymbiontS...]

HolobiontMetaPop = Container(MainEnvir([
    SubEnvir([ParamSubEnvir(
                           [SimParam.SubEnvir.Volume_subenvironments[sub]] # values
                           ,:ParamSubEnvir # Iam
                           ,SimParam.SubEnvir.Flow_sub_environements[sub] # Flow
                           ,SimParam.SubEnvir.ConcArriving[sub] # ConcArriving
                           ,[0.0 for i in 1:NallelesAggressiveness] # CurrentAgress
                           ,tuple([(Posi=Vector{Int}(),ID=Vector{Int}()) for sp in 1:SimParam.Nsp]...) # IndexNodesBySp
                           ,[2,sub] # PosiEnvir
                           )
             ,acellgenoS[  filter(i ->  SimParam.SubEnvir.PosiSp.ByEnvir[sub][i], eachindex(acellgenoS))...   ]
             ]
            ,[SimParam.SubEnvir.ConcArriving[sub]...].*SimParam.SubEnvir.Volume_subenvironments[sub]
            )
    for sub in 1:(SimParam.SubEnvir.NsubenvironmentsBySpecies*SimParam.Nsp)]
                                      )
                            )


for (i,sub) in enumerate(HolobiontMetaPop.nodes[2].nodes)
    PosiSp = HolobiontMetaPop.nodes[1].SimParam.SubEnvir.PosiSp.ByEnvir[i]
    push!(sub.nodes[1].IndexNodesBySp[findfirst(i ->i, PosiSp)].Posi,     2 )
    push!(sub.nodes[1].IndexNodesBySp[findfirst(i ->i, PosiSp)].ID  ,Int(1))
end

