library(MSLibrarian)

# create calibration library
diaFolder = "D:/O2TM_proteomics/1st_group"
diaFiles = list.files(diaFolder, pattern = ".raw$", full.names = T)
test_file = diaFiles[1]

projectFolder = "D:/O2TM_proteomics/Human_opt_spec_lib_MSLibrarian/"
fasta = "D:/O2TM_proteomics/Human_lib_with_DECM/human_uniprot_includeIsoform_reviewed_2022.11.10.fasta"
searchEngine = "msfragger"
irt = "biognosys_irt"

create.calibration.lib(projectFolder = projectFolder,
                        fasta = fasta,
                        diaFiles = test_file,
                        searchEngine = searchEngine,
                        irt = irt,
                        msfragger = 'C:/MSFragger-3.6/MSFragger-3.6.jar',
                        fragParams = 'C:/Fragparams/fragger_closed_mslibrarian_default.params', 
                        diaUmpireParams = 'C:/Fragparams/diaumpire_se_mslibrarian_default.params',
                        spectrastParams = 'C:/Fragparams/spectrast_create_mslibrarian_default.params',
                        msConvert = 'C:/Program Files/ProteoWizard/ProteoWizard 3.0.22317.1e024d4/msconvert.exe',
                        threads = 1)

predictionDb = "D:/O2TM_proteomics/human_prosit_intensity_hcd_2020_prosit_irt_2019.sqlite"
rt = "iRT"

process.calibration.lib(projectFolder = projectFolder,
                        predictionDb = predictionDb,
                        rt = rt)

system2('C:/Program Files/OpenMS-2.6.0/bin/DecoyDatabase.exe', args = c("-in",
                          fasta,
                          "-out",
                          '2.fasta'))
