import cobra.io
from cobra.io import read_sbml_model
from cobra import Reaction, Metabolite

def get_base_model():
    
    model = read_sbml_model('iNF517.xml')
    medium = model.medium
    original_medium = model.medium
    original_objective = model.objective
    

    # first new reaction, completion
    # https://www.kegg.jp/entry/R08165
    # http://bigg.ucsd.edu/universal/reactions/SEPHCHCS
    # Isochorismate + 2-Oxoglutarate <=> 2-Succinyl-5-enolpyruvyl-6-hydroxy-3-cyclohexene-1-carboxylate + CO2
    # akg_c + h_c + ichor_c ⇌ 2sephchc_c + co2_c, menD 

    model.remove_reactions(['SEPHCHCS'])
    new_reaction1 = Reaction('SEPHCHCS') # the enzyme / reaction name from BIGG
    new_reaction1.name = '2-succinyl-5-enolpyruvyl-6-hydroxy-3-cyclohexene-1-carboxylate synthase'
    # new_reaction1.gene_reaction_rule = '( )'
    sephchc = Metabolite(id='2sephchc_c', compartment='c', name='2-succinyl-5-enolpyruvyl-6-hydroxy-3-cyclohexene-1-carboxylate')
    new_reaction1.gene_reaction_rule = '( LLMG_1829 )'

    new_reaction1.add_metabolites({model.metabolites.get_by_id('akg_c'): -1,
                                   model.metabolites.get_by_id('h_c'): -1,
                                   model.metabolites.get_by_id('ichor_c'): -1,
                                   sephchc: 1,
                                   model.metabolites.get_by_id('co2_c'): 1
                                  })

    model.add_reactions([new_reaction1])

    # second new reaction, completion
    # https://www.kegg.jp/entry/R08166
    # http://bigg.ucsd.edu/universal/reactions/SHCHCS3
    # 2-Succinyl-5-enolpyruvyl-6-hydroxy-3-cyclohexene-1-carboxylate <=> 
    # (1R,6R)-6-Hydroxy-2-succinylcyclohexa-2,4-diene-1-carboxylate + Pyruvate
    # 2sephchc_c ⇌ 2shchc_c + pyr_c, menH

    model.remove_reactions(['SHCHCS3'])
    new_reaction2 = Reaction('SHCHCS3') # the enzyme / reaction name from BIGG
    new_reaction2.name = '2-succinyl-6-hydroxy-2,4-cyclohexadiene-1-carboxylate synthase'
    new_reaction2.gene_reaction_rule = '( LLMG_1830 )'

    new_reaction2.add_metabolites({model.metabolites.get_by_id('2sephchc_c'): -1,
                                   model.metabolites.get_by_id('2shchc_c'): 1,
                                   model.metabolites.get_by_id('pyr_c'): 1
                                  })

    model.add_reactions([new_reaction2])

    # third new reaction, completion
    # https://www.kegg.jp/entry/R07263
    # http://bigg.ucsd.edu/universal/reactions/SBZCOADH_x
    # 2-Succinylbenzoyl-CoA <=> 1,4-Dihydroxy-2-naphthoyl-CoA + H2O
    # h_c + sbzcoa_c ⇌ h2o_c + dhncoa_c, menB

    model.remove_reactions(['DHNCOAS'])
    new_reaction3 = Reaction('DHNCOAS') # the enzyme / reaction name from BIGG
    new_reaction3.name = '4-(2-carboxyphenyl)-4-oxobutanoyl-CoA dehydratase'
    new_reaction3.gene_reaction_rule = '( LLMG_1831 )'

    dhncoa = Metabolite(id='14dhncoa_c', compartment='c', name='1,4-Dihydroxy-2-naphthoyl-CoA')

    new_reaction3.add_metabolites({model.metabolites.get_by_id('sbzcoa_c'): -1,
                                   model.metabolites.get_by_id('h_c'): -1,
                                   dhncoa: 1,
                                   model.metabolites.get_by_id('h2o_c'): 1
                                  })

    model.add_reactions([new_reaction3])

    # fourth new reaction, completion
    # https://www.kegg.jp/entry/R07262
    # http://bigg.ucsd.edu/universal/reactions/NPHS_x
    # 1,4-Dihydroxy-2-naphthoyl-CoA + H2O <=> 1,4-Dihydroxy-2-naphthoate + CoA
    # h2o_c + 14dhncoa_c ⇌ coa_c + h_c + dhna_c, FabZ

    model.remove_reactions(['NPHS_c'])
    new_reaction4 = Reaction('NPHS_c') # the enzyme / reaction name from BIGG
    new_reaction4.name = '1,4-dihydroxy-2-naphthoyl-CoA hydrolase'
    new_reaction4.gene_reaction_rule = '( LLMG_1834 )'

    new_reaction4.add_metabolites({model.metabolites.get_by_id('h2o_c'): -1,
                                   model.metabolites.get_by_id('14dhncoa_c'): -1,
                                   model.metabolites.get_by_id('coa_c'): 1,
                                   model.metabolites.get_by_id('h_c'): 1,
                                   model.metabolites.get_by_id('dhna_c'): 1
                                  })

    model.add_reactions([new_reaction4])

    # fifth reaction, engineering or completion - ? In KEGG, there's just a reaciton, without the number (7). 
    # adding this reaction to obtain 2dmmq7_c from somewhere
    # not found in L. lactis though
    # https://www.kegg.jp/entry/R10757
    # http://bigg.ucsd.edu/universal/reactions/DHNAOT7
    # dhna_c + nad_c + hepdp_c ⇌ co2_c + nadh_c + ppi_c + 2dmmq7_c, menA

    model.remove_reactions(['DHNAOT7'])
    new_reaction5 = Reaction('DHNAOT7') 
    new_reaction5.name = '1,4-dihydroxy-2-naphthoate octaprenyltransferase'
    new_reaction5.gene_reaction_rule = '( LLMG_0197 )'

    dmmq7 = Metabolite(id='2dmmq7_c', compartment='c', name='2-Demethylmenaquinol-7')

    new_reaction5.add_metabolites({model.metabolites.get_by_id('dhna_c'): -1,
                                   model.metabolites.get_by_id('nad_c'): -1,
                                   model.metabolites.get_by_id('hepdp_c'): -1,
                                   model.metabolites.get_by_id('co2_c'): 1,
                                   model.metabolites.get_by_id('nadh_c'): 1,
                                   model.metabolites.get_by_id('ppi_c'): 1,
                                   dmmq7: 1
                                  })

    model.add_reactions([new_reaction5])

    # sixth reaction, completion
    # https://www.kegg.jp/entry/R09736
    # http://bigg.ucsd.edu/universal/reactions/AMMQT7
    # amet_c + nadph_c + 2dmmq7_c ⇌ ahcys_c + nadp_c + mql7_c, UbiE, MenG 

    model.remove_reactions(['AMMQT7'])
    new_reaction6 = Reaction('AMMQT7') 
    new_reaction6.name = 'Demethylmenaquinol methyltransferase'
    new_reaction6.gene_reaction_rule = '( LLMG_0753 )'


    new_reaction6.add_metabolites({model.metabolites.get_by_id('amet_c'): -1,
                                   model.metabolites.get_by_id('nadph_c'): -1,
                                   model.metabolites.get_by_id('2dmmq7_c'): -1,
                                   model.metabolites.get_by_id('ahcys_c'): 1,
                                   model.metabolites.get_by_id('nadp_c'): 1,
                                   model.metabolites.get_by_id('mql7_c'): 1
                                  })

    model.add_reactions([new_reaction6])

    # seventh reaction, engineering
    # https://www.kegg.jp/entry/R09991
    # http://bigg.ucsd.edu/universal/reactions/FADMQOR - not for mql7, but for mql8
    # fadh2_c + mqn7_c ⇌ fad_c + mql7_c

    new_reaction7 = Reaction('FADMQOR')
    new_reaction7.name = 'Menaquinone:FAD oxidoreductase' 

    new_reaction7.add_metabolites({model.metabolites.get_by_id('mql7_c'): -1,
                                   model.metabolites.get_by_id('fad_c'): -1,
                                   model.metabolites.get_by_id('mqn7_c'): 1,
                                   model.metabolites.get_by_id('fadh2_c'): 1
                                  })

    model.add_reactions([new_reaction7])

    # eighth reaction, engineering
    # http://bigg.ucsd.edu/universal/reactions/POX3
    # h2o_c + pyr_c + mqn7_c ⇌ ac_c + co2_c + mql7_c

    new_reaction8 = Reaction('POX3') 
    new_reaction8.name = 'Pyruvate:menaquinone oxidoreductase'

    new_reaction8.add_metabolites({model.metabolites.get_by_id('ac_c'): -1,
                                   model.metabolites.get_by_id('co2_c'): -1,
                                   model.metabolites.get_by_id('mql7_c'): -1,
                                   model.metabolites.get_by_id('h2o_c'): 1,
                                   model.metabolites.get_by_id('pyr_c'): 1,
                                   model.metabolites.get_by_id('mqn7_c'): 1
                                  })

    model.add_reactions([new_reaction8])

    # ninth reaction, engineering
    # http://bigg.ucsd.edu/universal/reactions/MQNS
    # https://www.uniprot.org/uniprotkb/Q9CF18/entry - enzyme in Lactococcus lactis IL1403
    # amet_c + 2dmmq7_c ⇌ ahcys_c + h_c + mqn7_c

    new_reaction9 = Reaction('MQNS')
    new_reaction9.name = '1,4-dihydroxy-2-naphthoate octaprenyltransferase'

    new_reaction9.add_metabolites({model.metabolites.get_by_id('amet_c'): -1,
                                   model.metabolites.get_by_id('2dmmq7_c'): -1,
                                   model.metabolites.get_by_id('ahcys_c'): 1,
                                   model.metabolites.get_by_id('h_c'): 1,
                                   model.metabolites.get_by_id('mqn7_c'): 1
                                  })

    model.add_reactions([new_reaction9])

    reaction_name = "MQNS"
    reaction = model.reactions.get_by_id(reaction_name)
    new_lower_bound = 0.0
    reaction.lower_bound = new_lower_bound
    model.reactions.MQNS.bounds
    model.add_boundary(model.metabolites.mqn7_c, type='demand')
    model.objective = model.reactions.DM_mqn7_c
    
    return model

def update_medium(model):

    # change medium
    original_medium = model.medium
    original_medium['EX_tyr__L_e'] = 0
    original_medium['EX_arg__L_e'] = 1
    original_medium['EX_glu__L_e'] = 16.05
    model.medium = original_medium

    return model

def add_sucrose(model):
    
    original_medium['EX_sucr_e'] = 284
    
    return model
    

