# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 08:11:56 2017

@author: phosp
"""

# the extraction function should only get platyhelminth orthologs, which requires this function
def is_platyhelminth(species):
    platys = ['clonorchis_sinensis_prjda72781',
              'diphyllobothrium_latum_prjeb1206',
              'echinococcus_canadensis_prjeb8992',
              'echinococcus_granulosus_prjeb121',
              'echinococcus_granulosus_prjna182977',
              'echinococcus_multilocularis_prjeb122',
              'echinostoma_caproni_prjeb1207',
              'fasciola_hepatica_prjeb6687',
              'fasciola_hepatica_prjna179522',
              'gyrodactylus_salaris_prjna244375',
              'hydatigera_taeniaeformis_prjeb534',
              'hymenolepis_diminuta_prjeb507',
              'hymenolepis_microstoma_prjeb124',
              'hymenolepis_nana_prjeb508',
              'macrostomum_lignano_prjna284736',
              'mesocestoides_corti_prjeb510',
              'opisthorchis_viverrini_prjna222628',
              'protopolystoma_xenopodis_prjeb1201',
              'schistocephalus_solidus_prjeb527',
              'schistosoma_curassoni_prjeb519',
              'schistosoma_haematobium_prjna78265',
              'schistosoma_japonicum_prjea34885',
              'schistosoma_mansoni_prjea36577',
              'schistosoma_margrebowiei_prjeb522',
              'schistosoma_mattheei_prjeb523',
              'schistosoma_rodhaini_prjeb526',
              'schmidtea_mediterranea_prjna12585',
              'spirometra_erinaceieuropaei_prjeb1202',
              'taenia_asiatica_prjeb532',
              'taenia_asiatica_prjna299871',
              'taenia_saginata_prjna71493',
              'taenia_solium_prjna170813',
              'trichobilharzia_regenti_prjeb4662']
    logic = species in platys    
    return(logic)
