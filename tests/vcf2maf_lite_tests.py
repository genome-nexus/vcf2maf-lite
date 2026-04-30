#!/usr/bin/env python3

import unittest
import tempfile
import os
from vcf2maf_lite.vcf2maf_lite import extract_vcf_data_from_file
from vcf2maf_lite.vcf2maf_lite import write_standardized_mutation_file

class vcf2maf_lite_DataTests(unittest.TestCase):

    def setUp(self):
        self.temp_files = []

    def tearDown(self):
        for temp_file in self.temp_files:
            if os.path.exists(temp_file):
                os.remove(temp_file)

    def test_extract_vcf_data_from_file_no_samples(self):
        _, vcf = tempfile.mkstemp()
        self.temp_files.append(vcf)
        with open(vcf, 'w') as f:
           f.write(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n"
            "20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\n"
           )
        with self.assertRaises(Exception) as exc:
            extract_vcf_data_from_file(vcf, 'center name 1', 'sequence source 1')
        self.assertEqual("No sample column found", str(exc.exception))

    def test_extract_vcf_data_from_file_1_sample(self):
        _, vcf = tempfile.mkstemp()
        self.temp_files.append(vcf)
        with open(vcf, 'w') as f:
           f.write(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
            "20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\n"
           ) 
        maf_data = extract_vcf_data_from_file(vcf, 'center name 1', 'sequence source 1')
        self.assertEqual(1, len(maf_data))
        maf_row = maf_data[0]
        self.assertEqual('S1', maf_row['Tumor_Sample_Barcode']) 
        self.assertEqual('NORMAL', maf_row['Matched_Norm_Sample_Barcode']) 

    def test_extract_vcf_data_from_file_2_samples(self):
        _, vcf = tempfile.mkstemp()
        self.temp_files.append(vcf)
        with open(vcf, 'w') as f:
           f.write(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
            "20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51\n"
           ) 
        maf_data = extract_vcf_data_from_file(vcf, 'center name 1', 'sequence source 1')
        self.assertEqual(1, len(maf_data)) 
        maf_row = maf_data[0]
        self.assertEqual('S1', maf_row['Tumor_Sample_Barcode']) 
        self.assertEqual('S2', maf_row['Matched_Norm_Sample_Barcode']) 

    def test_extract_vcf_data_from_file_tumor(self):
        _, vcf = tempfile.mkstemp()
        self.temp_files.append(vcf)
        with open(vcf, 'w') as f:
           f.write(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\n"
            "20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\n"
           )
        maf_data = extract_vcf_data_from_file(vcf, 'center name 1', 'sequence source 1')
        self.assertEqual(1, len(maf_data)) 
        maf_row = maf_data[0]
        self.assertTrue(maf_row['Tumor_Sample_Barcode']) 
        self.assertNotEqual('S1', maf_row['Tumor_Sample_Barcode']) 
        self.assertEqual('NORMAL', maf_row['Matched_Norm_Sample_Barcode']) 

    def test_extract_vcf_data_from_file_3_samples(self):
        _, vcf = tempfile.mkstemp()
        self.temp_files.append(vcf)
        with open(vcf, 'w') as f:
           f.write(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n"
            "20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51\t1|0:48:8:51,51\n"
           )
        with self.assertRaises(Exception) as exc:
            extract_vcf_data_from_file(vcf, 'center name 1', 'sequence source 1')
        self.assertEqual("Expected max 2 sample columns for tumor and normal sample. But found 3 columns.", str(exc.exception))

    def test_extract_vcf_data_from_file_tumor_normal_swap(self):
        # this test is just showing that the order of the columns does not matter...
        _, vcf = tempfile.mkstemp()
        self.temp_files.append(vcf)
        tumor_sample_id = os.path.basename(vcf)
        with open(vcf, 'w') as f:
           f.write(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR\n"
            "20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51\n"
           )
        maf_data = extract_vcf_data_from_file(vcf, 'center name 1', 'sequence source 1')
        self.assertEqual(1, len(maf_data))
        maf_row = maf_data[0]
        self.assertEqual(tumor_sample_id, maf_row['Tumor_Sample_Barcode'])
        self.assertEqual('NORMAL', maf_row['Matched_Norm_Sample_Barcode'])

    def test_extract_vcf_data_from_file_1_normal_sample(self):
        _, vcf = tempfile.mkstemp()
        self.temp_files.append(vcf)
        with open(vcf, 'w') as f:
           f.write(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\n"
            "20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\n"
           )
        with self.assertRaises(Exception) as exc:
            extract_vcf_data_from_file(vcf, 'center name 1', 'sequence source 1')
        self.assertEqual("There is only one sample column and it has NORMAL label. No tumor sample column present.", str(exc.exception))

    def test_extract_vcf_data_from_file_2_samples_specified_in_header(self):
        _, vcf = tempfile.mkstemp()
        self.temp_files.append(vcf)
        with open(vcf, 'w') as f:
           f.write(
            '##normal_sample=S1\n'
            '##tumor_sample=S2\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
            "20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51\n"
           )
        maf_data = extract_vcf_data_from_file(vcf, 'center name 1', 'sequence source 1')
        self.assertEqual(1, len(maf_data))
        maf_row = maf_data[0]
        self.assertEqual('S2', maf_row['Tumor_Sample_Barcode'])
        self.assertEqual('S1', maf_row['Matched_Norm_Sample_Barcode'])

    def test_extract_vcf_data_from_file_2_samples_specified_in_header_reversed_order(self):
        # same as above, but with column order reversed in CHROM header line...should still give the same result, just
        # to show that it is really reading from header and not from column order:
        _, vcf = tempfile.mkstemp()
        self.temp_files.append(vcf)
        with open(vcf, 'w') as f:
           f.write(
            '##normal_sample=S1\n'
            '##tumor_sample=S2\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS2\tS1\n"
            "20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51\n"
           )
        maf_data = extract_vcf_data_from_file(vcf, 'center name 1', 'sequence source 1')
        self.assertEqual(1, len(maf_data))
        maf_row = maf_data[0]
        self.assertEqual('S2', maf_row['Tumor_Sample_Barcode'])
        self.assertEqual('S1', maf_row['Matched_Norm_Sample_Barcode'])

    def test_extract_vcf_data_from_file_normal_sample_refers_to_non_existing_column(self):
        _, vcf = tempfile.mkstemp()
        self.temp_files.append(vcf)
        with open(vcf, 'w') as f:
           f.write(
            '##normal_sample=S1\n'
            '##tumor_sample=S2\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS2\n"
            "20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\n"
           )
        with self.assertRaises(Exception) as exc:
            extract_vcf_data_from_file(vcf, 'center name 1', 'sequence source 1')
        self.assertEqual("There is normal_sample=S1 in the header, but no respective column found.", str(exc.exception))

    def test_extract_vcf_data_from_file_tumor_sample_refers_to_non_existing_column(self):
        _, vcf = tempfile.mkstemp()
        self.temp_files.append(vcf)
        with open(vcf, 'w') as f:
           f.write(
            '##normal_sample=S1\n'
            '##tumor_sample=S2\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
            "20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\n"
           )
        with self.assertRaises(Exception) as exc:
            extract_vcf_data_from_file(vcf, 'center name 1', 'sequence source 1')
        self.assertEqual("There is tumor_sample=S2 in the header, but no respective column found.", str(exc.exception))

    def test_extract_vcf_data_from_file_header_samples_have_precedence_over_the_rest_of_strategies(self):
        _, vcf = tempfile.mkstemp()
        self.temp_files.append(vcf)
        with open(vcf, 'w') as f:
           f.write(
            '##normal_sample=TUMOR\n'
            '##tumor_sample=NORMAL\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\tNORMAL\n"
            "20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51\n"
           )
        maf_data = extract_vcf_data_from_file(vcf, 'center name 1', 'sequence source 1')
        self.assertEqual(1, len(maf_data))
        maf_row = maf_data[0]
        self.assertEqual('NORMAL', maf_row['Tumor_Sample_Barcode'])
        self.assertEqual('TUMOR', maf_row['Matched_Norm_Sample_Barcode'])

    def test_extract_vcf_data_from_file_multiple_equals_in_header(self):
        _, vcf = tempfile.mkstemp()
        self.temp_files.append(vcf)
        with open(vcf, 'w') as f:
           f.write(
            '##tumor_sample=A=B\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA=B\n"
            "20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\n"
           )
        with self.assertRaises(Exception) as exc:
            extract_vcf_data_from_file(vcf, 'center name 1', 'sequence source 1')
        self.assertIn("The tumor_sample and normal_sample are expected together", str(exc.exception))
        
    def test_extract_vcf_data_from_file_invalid_data_file(self):
        _, vcf = tempfile.mkstemp()
        self.temp_files.append(vcf)
        with open(vcf, 'w') as f:
            f.write(
            '##normal_sample=TUMOR\n'
            '##tumor_sample=NORMAL\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\tNORMAL\n"
            )
        maf_data = extract_vcf_data_from_file(vcf, 'center name 1', 'sequence source 1')
        self.assertEqual([], maf_data)
    
    def test_extract_vcf_data_from_file_malformed_records(self):
        _, vcf = tempfile.mkstemp()
        self.temp_files.append(vcf)
        with open(vcf, 'w') as f:
            f.write(
            '##normal_sample=TUMOR\n'
            '##tumor_sample=NORMAL\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\tNORMAL\n"
            "20\t14370\trs6054257\tG\tA\t29\tPASS\t\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51\n"
            )
        maf_data = extract_vcf_data_from_file(vcf, 'center name 1', 'sequence source 1')
        self.assertEqual(None, maf_data)
    
    def test_extract_vcf_data_from_file_varscan_allele_counts(self):
        _, vcf = tempfile.mkstemp()
        self.temp_files.append(vcf)
        with open(vcf, 'w') as f:
            f.write(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
            "1\t10105\t.\tA\tC\t7\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:AD:RD:FT:GQ:GL\t1/0:30:15:PASS:13:-1.58548,-0.0193946,-153.729\n"
            )
        maf_data = extract_vcf_data_from_file(vcf, 'center name 1', 'sequence source 1')
        self.assertEqual(1, len(maf_data))
        maf_row = maf_data[0]
        self.assertEqual('15', maf_row['t_ref_count'])
        self.assertEqual('30', maf_row['t_alt_count'])
        self.assertEqual('45', maf_row['t_depth'])
        
    def test_extract_vcf_data_from_file_strelka_snp_allele_counts(self):
        _, vcf = tempfile.mkstemp()
        self.temp_files.append(vcf)
        with open(vcf, 'w') as f:
            f.write(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\n"
            "1\t10105\t.\tA\tN\t7\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:AU:CU:GU:TU:FT\t1/0:20:15:10:10:PASS\n"
            )
        maf_data = extract_vcf_data_from_file(vcf, 'center name 1', 'sequence source 1')
        self.assertEqual(1, len(maf_data))
        maf_row = maf_data[0]
        self.assertEqual('20', maf_row['t_ref_count'])
        self.assertEqual('15', maf_row['t_alt_count'])
        self.assertEqual('35', maf_row['t_depth'])
        
    def test_extract_vcf_data_from_file_bcftools_allele_counts(self):
        _, vcf = tempfile.mkstemp()
        self.temp_files.append(vcf)
        with open(vcf, 'w') as f:
            f.write(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\n"
            "1\t10105\t.\tA\tT\t7\tPASS\tNS=3;DP=50;AF=0.5;DB;H2\tGT:DV:DP:FT\t1/0:30:50:PASS\n"
            )
        maf_data = extract_vcf_data_from_file(vcf, 'center name 1', 'sequence source 1')
        self.assertEqual(1, len(maf_data))
        maf_row = maf_data[0]
        self.assertEqual('30', maf_row['t_alt_count'])
        self.assertEqual('50', maf_row['t_depth'])
        
    def test_extract_vcf_data_from_file_resolve_alleles_positions_variant_class(self):
        _, vcf = tempfile.mkstemp()
        self.temp_files.append(vcf)
        with open(vcf, 'w') as f:
            f.write(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR\n"
            "1\t45796859\t.\tTCATGGCGGTGG\tT\t.\tPASS\tCONTQ=93;DP=1307;ECNT=1\tGT:AD:AF:DP\t0/0:343,1:0.005711:344\t0/1:920,5:0.006427:925\n"
            "1\t10496754\t.\tA\tACC\t.\tPASS\tDP=104;ECNT=1\tGT:AD:AF:DP\t1/0:34,23:0.005711:344\t0/1:920,5:0.006427:925\n"
            )
        maf_data = extract_vcf_data_from_file(vcf, 'center name 1', 'sequence source 1')
        self.assertEqual(2, len(maf_data))
        maf_row = maf_data[0]
        self.assertEqual('Frame_Shift_Del', maf_row['Variant_Classification'])
        self.assertEqual('DEL', maf_row['Variant_Type'])
        self.assertEqual('CATGGCGGTGG', maf_row['Reference_Allele'])
        self.assertEqual('-', maf_row['Tumor_Seq_Allele2'])
        self.assertEqual('45796860', maf_row['Start_Position'])
        self.assertEqual('45796870', maf_row['End_Position'])

        # Verify insertion coordinates: VCF POS should NOT be incremented for insertions
        maf_row_ins = maf_data[1]
        self.assertEqual('Frame_Shift_Ins', maf_row_ins['Variant_Classification'])
        self.assertEqual('INS', maf_row_ins['Variant_Type'])
        self.assertEqual('-', maf_row_ins['Reference_Allele'])
        self.assertEqual('CC', maf_row_ins['Tumor_Seq_Allele2'])
        self.assertEqual('10496754', maf_row_ins['Start_Position'])
        self.assertEqual('10496755', maf_row_ins['End_Position'])

    def test_insertion_coordinates_match_vcf2maf_pl(self):
        """Test for GitHub issue #6: insertion coordinates should match vcf2maf.pl output."""
        _, vcf = tempfile.mkstemp()
        self.temp_files.append(vcf)
        with open(vcf, 'w') as f:
            f.write(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\n"
            "1\t116579677\t.\tC\tCCTA\t.\tPASS\tDP=100\tGT:AD:DP\t0/1:50,50:100\n"
            )
        maf_data = extract_vcf_data_from_file(vcf, 'center name 1', 'sequence source 1')
        self.assertEqual(1, len(maf_data))
        maf_row = maf_data[0]
        self.assertEqual('-', maf_row['Reference_Allele'])
        self.assertEqual('CTA', maf_row['Tumor_Seq_Allele2'])
        self.assertEqual('INS', maf_row['Variant_Type'])
        # Per vcf2maf.pl expected output: start=116579677, end=116579678
        self.assertEqual('116579677', maf_row['Start_Position'])
        self.assertEqual('116579678', maf_row['End_Position'])
        

    def test_extract_vcf_data_from_file_na_in_output_file(self):
        _, vcf = tempfile.mkstemp()
        self.temp_files.append(vcf)
        with open(vcf, 'w') as f:
            f.write(
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR\n"
                "1\t45796859\t.\tTCATGGCGGTGG\tT\t.\tPASS\tMuTect2;Custom_filters=strand_bias;RepeatMasker=Simple_repeat;DP=1191\tGT:AD:AF:DP:alt_count_raw:ref_count_raw:depth_raw\t0/0:343,1:0.005711:344:0:424:425\t0/1:920,5:0.006427:925:0:1184:1191\n"
                "2\t114020823\t.\tG\tGT\t.\tPASS\tStrelka2;Custom_filters=strand_bias;RepeatMasker=Low_complexity;PoN=53;DP=53\tGT:AD:AF:DP:alt_count_raw:ref_count_raw:depth_raw\t1/0:34,23:0.005711:344:1:12:13\t0/1:920,5:0.006427:925:7:46:53\n"
                "6\t32521751\t.\tA\tT\t.\tPASS\tMuTect2;Strelka2;Strelka2FILTER;DP=677\tGT:AD:AF:DP:alt_count_raw:ref_count_raw:depth_raw\t0/0:343,1:0.005711:344:0:159:160\t1/0:34,23:0.005711:344:20:657:677\n"
            )
        maf_data = extract_vcf_data_from_file(vcf, 'center name 1', 'sequence source 1', 'TUMOR', 'NORMAL', ['MuTect2','Strelka2','Custom_filters','Strelka2FILTER','RepeatMasker','PoN'], ['alt_count_raw','ref_count_raw','depth_raw'])
        self.assertEqual(3, len(maf_data))
        maf_row = maf_data[0]
        self.assertEqual('1', maf_row['MuTect2'])
        self.assertEqual('NA', maf_row['Strelka2'])
        self.assertEqual('1184', maf_row['t_ref_count_raw'])
        self.assertEqual('425', maf_row['n_depth_raw'])

        with tempfile.NamedTemporaryFile(delete=False) as temp_file:
            self.temp_files.append(temp_file.name)
            write_standardized_mutation_file(maf_data, temp_file.name)

        with open(self.temp_files[-1], 'r') as f:
            maf_content = f.read()

        self.assertIn('\tNA\t', maf_content)

if __name__=='__main__':
    unittest.main()
