library(reticulate)
use_python("/usr/bin/python3")
library(SigProfilerExtractorR)

samples = "Our_Zou_Kucab_SNVs.txt"
sigprofilerextractor(input_type = "matrix",
                    output = "Signatures/All_samples",
                    input_data = samples, reference_genome="GRCh38", opportunity_genome="GRCh38",
                    minimum_signatures=1, maximum_signatures=15, nmf_replicates=100, cpu=8, gpu=F, batch_size=1)

samples = "Our_Zou_Kucab_SNVs_sub1.txt"
sigprofilerextractor(input_type = "matrix",
                   output = "Signatures/Sub1",
                    input_data = samples, reference_genome="GRCh38", opportunity_genome="GRCh38",
                    minimum_signatures=1, maximum_signatures=10, nmf_replicates=100, cpu=8, gpu=F, batch_size=1)

samples = "Our_Zou_Kucab_SNVs_sub2.txt"
sigprofilerextractor(input_type = "matrix",
                    output = "Signatures/Sub2",
                    input_data = samples, reference_genome="GRCh38", opportunity_genome="GRCh38",
                    minimum_signatures=1, maximum_signatures=10, nmf_replicates=100, cpu=8, gpu=F, batch_size=1)

samples = "Our_Zou_Kucab_SNVs_sub3.txt"
sigprofilerextractor(input_type = "matrix",
                    output = "Signatures/Sub3",
                    input_data = samples, reference_genome="GRCh38", opportunity_genome="GRCh38",
                    minimum_signatures=1, maximum_signatures=10, nmf_replicates=100, cpu=8, gpu=F, batch_size=1)

samples = "Our_Zou_Kucab_indels.txt"
sigprofilerextractor(input_type = "matrix", 
                     output = "Signatures/All_samples",
                     input_data = samples, reference_genome="GRCh38", opportunity_genome="GRCh38", 
                     minimum_signatures=1, maximum_signatures=10, nmf_replicates=100, cpu=8, gpu=F, batch_size=1)
