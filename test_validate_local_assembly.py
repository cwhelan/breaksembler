import unittest
import tempfile
import validate_local_assembly
import shutil
import pybedtools

class SVsFromSAMTest(unittest.TestCase):
    def test_deletion_in_alignment(self):
        # need to create temp files here until I can figure out how to mock out pysam
        temp_dir = tempfile.mkdtemp()
        
        try :
            sam_records = \
"""@SQ	SN:2	LN:242951149
NODE_1_length_2155_cov_7.27857_ID_3018	16	2	168788932	60	1050M76D1006M1I98M	*	0	0	AAAAGAATAGGAAATTCTTAAACTAAAAAAACAAACAACAAAAAAAAATCGCTTTCTTAGAGCACTGTTATGCTACAAGATAACTGATGAAGAGGCTCTATACATACATATGTACTGATAATAATGCTGATAATAATGCTTATAAAGATCCATCAGCATTCAAATAAAAAAGGATAAACAAGGGGACACAGTCAAAGGGAATATGGTTCAACTGTTTATTTAAATCTGGTAATCAGATAAACCATCAATCTATAACAAAATAACAGTAGGAAAATTCATTTCATTATTGCATGAAATGAAAACCAAAAGAAAGCGCATTTACAAATGCACACTAAGCCAGGCATTAGGTGATTTTCAACTCCTCACCAAAGCAATCAAATAGAGATTAACATACCTGTCAGCTGAAACCAGGCTGTCAAAATGAAAGCCTTAAATTACACATGTCCTTAACAAAAAGGATTTAAGAGATTTTTCTCCTAACTAATAAGGAGCACTGCATAGTTTTTAATAAAACAATAATGTTTAATTTTGCTTATGTGCTAAAAAGGCATTTTTAAGGAGATTAAAACGCAGATTCAAGACCACCAGAAGAAACTTTCCAACATGAAATCAGACCATTTGGTCATCTCTTCTGAGTACGTTGTTCATTTGTCTTTTACCCAGCAACAATTCCTGCAAAATGATAATCAAAACAGCCTTTTGGTAAACAGGAGCTTAAACACTGTCTCTCAAATCAGATTATATTTCAAGTGATCTGAAACTAAGACAGCACTTCATGAACATGGTTCTCATTTTAAAAAGTCACAGGGAAACTGCCCTTTGAGGTTCTGAGGGACACAGACACCAAAGTTTAGGACTCAGCAGTTCATGATCACTGACATATCCAACCCTCTAGCCATCGTTTTGAAAGTGAAGCAGAATGGGGTTCAGGAGGAAGAGAACAGGGCAAAAAGGGATACTTAATAGTGACCCTGAAGAGGAAAACAAAGCAACCTTGAGAAATCTTCCAGGGGATATTAACAAACTTGCATCTCTGAAACACTGGACCACCAAATACATCATAAAATTTGGGAGGATCAGATTCTCTTCTCAATATCCATGGAGAAGATACTGCTGATAAAACTATTGACATAGGAAAGTTAGTTCCTCAAAGAGTATCAGGTTCCAGGTATATTTTTCAGAATTTTATCACCATTACATAAACAGGATGAACTACTTCATAGTACTTTAAATCCAGCAGAAAAGGCAAGACCAGACCATATCTCAGGATTCGGGTATCATTATTTCAGAATTTAAGACAATTTATAACTGAAACTGGGCATAAGATTTCCCAATTTACAATCCATAAGGGAAAAGGGAGGTTCTTCCTAGGGAAGAAAGTTATTTATTGACTCAAAACAATAAACTGCTCTTTAAACACATTGTAAGTATATATTTTAACAAAGAAAGTACATTAGATCAGGATTACATTTTGCTGATTAGGGTGACATGTCGTCTCTGTCTCGTGACCAGTTACTTACCTAAGAGCAGTATTACATAAAACACTCACTAGTGGTGCCTGCCACACTTGTGATTCCTTCCCAGGCAGGAACTTCCGATGCTACACCCTTCTCCAAGGAGGGTGTCTGAGTTTAGGGCAGTCAGGTGTTCTATGAAGGGACTCACTAAAAACTAGATTCTCAATCTTTGGGATTAAGATTGGGAAATACAGACAGAATCCAGCAATTAGCAGCCATAGTTGCGGCTAACAGGGAAAGAAGAGAGAAGCTGGATGGGCCATAAGGAGCTATGAGTGAGCCAAAATTCTGACTAAGTAGAAACAACATGACAGAGAAGACATGCACAGATAGTAGGTAAGAAAGCCAGGTAAAAACAGCCAGAGAAGCAGGTGGGCATGAGCAAAGCAAACACATGCTGATTGGCCGGTTACACTAAGGTTAGAGGATCCAGAATAGTGAGGATCCAAAGTTTCCTGCTGCTGAGGTCTTATCCACTTCCCAGTCCATATGCAGACATTTCACTGCAATTCCCAAGACCGAGACCCTCAACTTTTCCAGGGATGCTAGAAAACCCAGGATTATGGATGAGATTCCTGTCCCCCTATTCCCACATCTACCTTTATAAAGAATGCCCCATTATTTGGGAACTTTTAGTG	*	NM:i:81	AS:i:2045	XS:i:0
"""
            sam_file_name = temp_dir + "/test.sam"
            sam_file = open(sam_file_name, "w")
            sam_file.write(sam_records)
            sam_file.close()
            print sam_file_name
            (samfile, records_by_name) = validate_local_assembly.read_sam_file(sam_file_name)
            
            deletion_intervals = validate_local_assembly.get_deletion_intervals(samfile, records_by_name)
            self.assertEqual(1, len(deletion_intervals))
            expected_deletion = pybedtools.Interval("2", 168789981, 168790057, "NODE_1_length_2155_cov_7.27857_ID_3018-1")
            print deletion_intervals[0]
            print expected_deletion
            self.assertEqual(expected_deletion, deletion_intervals[0])
        finally:
            shutil.rmtree(temp_dir)

    def test_deletion_between_alignments(self):
        # need to create temp files here until I can figure out how to mock out pysam
        temp_dir = tempfile.mkdtemp()
        
        try :
            sam_records = \
"""@SQ	SN:2	LN:242951149
NODE_1_length_2098_cov_6.69848_ID_3108	16	2	233375525	60	1105M993S	*	0	0	CATTCTAGGAATAAATAAATCCCACTTGGATTATGGTGTATAATCTTTTTAATATATTGCTGAATTTGCTTTGCTAGTATTTTGTTAAGAATTTTTGCATCAGTGTTCATAAAGGATATTGGTCTATATTTTTCTTTAGAGCCTTTGTCTGACTTTGATACCAGGGTAATTTTGGCTTCATAGAATTAGTTAGGAAGTGTAGGTAATCTTGCTTCCTTTTCAGTTTTTGGAAAAAGTTTGAGAAGGATTATGATTAGTTTTTCCTTAAATGTTTGGTGCGTTTTTCACCAGAGAAGCCATCAGGGCAGGGCTTTCCTTTGTCCAGAGATTCTTTTTTTTTTTTTTTAAAGACAGAGTCTAGCTCTGTCGGCCAGGCTGGAGTGCAGTGGCACAGTCTCAGCCCACTGCAACCGCTGTCTCCCGAGCTCAAGCAGTTCTCCTGCCTCAGCCTCCTGAGTAGCTGAGATTACAGGCGTGTGCCACCATGCCCGGCTCATTTTAGTAGAGATGGGGTTTCACCATGTTGGCCAGGCTGGTCTCGAACTCCTGACCTCAGATAATCTGCCTGCCTTGGCCTCCCAAACTGCTGGGATTACTGGCGTGAGCCACTGTGCCCGGCCTGTTCAGAGATTTTTGATTACTGATTCTAGATACTAGTCCTTTGATGGATATGTGGTTTGTACATATTTTTTGTCAGTCTTTATCTTTTCATTTCATTCTCTTAATAGAGTCTTTCACAGAGCAAAAAGGTGTAAATTTTGATTAATCCAATAGATCAATGTTTTTCTTTTATAGGTTGTTCTTTTTTGTATTAAGTTTAACAAGAACACTTTGCCTAGCTTTAGATGCCAATGGTTTTTTTTGTGTGTGTGTGGAAAAATTTTATAGTTTTTATAATTTACATTTAAGTTCATGGTCCATTTTAAATTAATTTTCATATAAGCTAGTTGTTTTTTGTTGTTGTTTTCTCAATGAGTGTCCGTTTGCTTCGGGACCAAGAATGTCTTTTCTCTCTTAAATTGTTTTCTTTCACCTTTGTCAAAATTTAGTTAAGCATATTTGTGTGGGTCTATTTCTTTGTGGGTCTGTTTCTATTTCTGCTTCTAATTTTGTCCCATAGATCACTGTCTCTTCCTTCACCAGTACTGTACAATCTTGATTACTGTAGCTGTATAATAAGTCTTGAAATTGGGTAGACTGATTCTTCCCAGTTTAATCTTCTTTTAAAAAAATTTTTTTTTAGCTTTGTATTTCCTTCGCCTTTCCATAAAAATTTTTAGGATAATTTTGTCCGTATTTACTAAAAATCTTGCTGAGATTCTGATAGGAATTATGTTAAACCTCTCTATCAGTTTGAGAAGACTTGATGTCTTTACCATGTTAGTCTTTTAGTTCATGAACTTGGTATGTCTATCTGTGTTGTTAGATCTTTGATTTCTTTCATCAGTGTTGCGTAGTTTTCAGCATGCATTCTGTACATGTAGATTTACATCTGAGCATCTCAGTTTTTTTTTTTTTTTTTCAGTGATAGTAAATGAGCTCATTTATTAGCTCAAGGAGTTCTTTTTGTTTTCTGTTTTTTTGGGTAGTTTCCTTGAGATTTTACATAGACAGTCATATTTATCCTTTTTTTGAGATGGAGTCTCACTCTGTCACCCAGGCTGGAGTGCAGTGGCGCGATCTCGGCTCCCTGCAACCTCTGCTTCCCAGGTCCAAGCGATTCTCCTGCCTCAGCCTCCCAAGTAGCTGGGACTACAGGCATGCGCCACCATGCCCAGCTAATTTTTGCATTTTTTAGTAGAGACAAGGTTTCACCATGTTGGCCAGGCTGGTCTCGAACTCCTGACCTCAGGTGATCTACCCACCTTGGCCTCCTAAAATGCTGGAATTAGAGGCGTGAACTGCTGTGCTCAGCCTGTATTTTTATTTTCATTCCGTTGTGTGTATTTTAATTTCCCTTGAGACTCCCTCTTTTACCCATGGATTATTTTGAAGTATATTTTTGGGTTTCTAAGTATTTGGAAGATTTCTGTTATTTTTCTTTTATTGATTTCTAGCTTGATTCTTTGTGGTCAAAGAAACATACTATATGTGATTT	*	NM:i:2	AS:i:1095	XS:i:43	SA:Z:2,233376921,-,1091S13M1D994M,60,2;
NODE_1_length_2098_cov_6.69848_ID_3108	2064	2	233376921	60	1091H13M1D994M	*	0	0	CTATTTCTGCTTCTAATTTTGTCCCATAGATCACTGTCTCTTCCTTCACCAGTACTGTACAATCTTGATTACTGTAGCTGTATAATAAGTCTTGAAATTGGGTAGACTGATTCTTCCCAGTTTAATCTTCTTTTAAAAAAATTTTTTTTTAGCTTTGTATTTCCTTCGCCTTTCCATAAAAATTTTTAGGATAATTTTGTCCGTATTTACTAAAAATCTTGCTGAGATTCTGATAGGAATTATGTTAAACCTCTCTATCAGTTTGAGAAGACTTGATGTCTTTACCATGTTAGTCTTTTAGTTCATGAACTTGGTATGTCTATCTGTGTTGTTAGATCTTTGATTTCTTTCATCAGTGTTGCGTAGTTTTCAGCATGCATTCTGTACATGTAGATTTACATCTGAGCATCTCAGTTTTTTTTTTTTTTTTTCAGTGATAGTAAATGAGCTCATTTATTAGCTCAAGGAGTTCTTTTTGTTTTCTGTTTTTTTGGGTAGTTTCCTTGAGATTTTACATAGACAGTCATATTTATCCTTTTTTTGAGATGGAGTCTCACTCTGTCACCCAGGCTGGAGTGCAGTGGCGCGATCTCGGCTCCCTGCAACCTCTGCTTCCCAGGTCCAAGCGATTCTCCTGCCTCAGCCTCCCAAGTAGCTGGGACTACAGGCATGCGCCACCATGCCCAGCTAATTTTTGCATTTTTTAGTAGAGACAAGGTTTCACCATGTTGGCCAGGCTGGTCTCGAACTCCTGACCTCAGGTGATCTACCCACCTTGGCCTCCTAAAATGCTGGAATTAGAGGCGTGAACTGCTGTGCTCAGCCTGTATTTTTATTTTCATTCCGTTGTGTGTATTTTAATTTCCCTTGAGACTCCCTCTTTTACCCATGGATTATTTTGAAGTATATTTTTGGGTTTCTAAGTATTTGGAAGATTTCTGTTATTTTTCTTTTATTGATTTCTAGCTTGATTCTTTGTGGTCAAAGAAACATACTATATGTGATTT	*	NM:i:2	AS:i:995	XS:i:0	SA:Z:2,233375525,-,1105M993S,60,2;
"""
            sam_file_name = temp_dir + "/test.sam"
            sam_file = open(sam_file_name, "w")
            sam_file.write(sam_records)
            sam_file.close()
            print sam_file_name
            (samfile, records_by_name) = validate_local_assembly.read_sam_file(sam_file_name)
            
            deletion_intervals = validate_local_assembly.get_deletion_intervals(samfile, records_by_name)
            print deletion_intervals
            self.assertEqual(2, len(deletion_intervals))
            expected_deletion1 = pybedtools.Interval("2", 233376629, 233376920, "NODE_1_length_2155_cov_7.27857_ID_3018-1")
            expected_deletion2 = pybedtools.Interval("2", 233376933, 233376934, "NODE_1_length_2155_cov_7.27857_ID_3018-2")
            print deletion_intervals[0]
            print expected_deletion1
            self.assertEqual(expected_deletion1, deletion_intervals[0])
            print deletion_intervals[1]
            print expected_deletion2
            self.assertEqual(expected_deletion2, deletion_intervals[1])
        finally:
            shutil.rmtree(temp_dir)

    def test_deletion_between_alignments(self):
        # need to create temp files here until I can figure out how to mock out pysam
        temp_dir = tempfile.mkdtemp()
        
        try :
            sam_records = \
"""@SQ	SN:2	LN:242951149
NODE_1_length_2098_cov_6.69848_ID_3108	16	2	233375525	60	1105M993S	*	0	0	CATTCTAGGAATAAATAAATCCCACTTGGATTATGGTGTATAATCTTTTTAATATATTGCTGAATTTGCTTTGCTAGTATTTTGTTAAGAATTTTTGCATCAGTGTTCATAAAGGATATTGGTCTATATTTTTCTTTAGAGCCTTTGTCTGACTTTGATACCAGGGTAATTTTGGCTTCATAGAATTAGTTAGGAAGTGTAGGTAATCTTGCTTCCTTTTCAGTTTTTGGAAAAAGTTTGAGAAGGATTATGATTAGTTTTTCCTTAAATGTTTGGTGCGTTTTTCACCAGAGAAGCCATCAGGGCAGGGCTTTCCTTTGTCCAGAGATTCTTTTTTTTTTTTTTTAAAGACAGAGTCTAGCTCTGTCGGCCAGGCTGGAGTGCAGTGGCACAGTCTCAGCCCACTGCAACCGCTGTCTCCCGAGCTCAAGCAGTTCTCCTGCCTCAGCCTCCTGAGTAGCTGAGATTACAGGCGTGTGCCACCATGCCCGGCTCATTTTAGTAGAGATGGGGTTTCACCATGTTGGCCAGGCTGGTCTCGAACTCCTGACCTCAGATAATCTGCCTGCCTTGGCCTCCCAAACTGCTGGGATTACTGGCGTGAGCCACTGTGCCCGGCCTGTTCAGAGATTTTTGATTACTGATTCTAGATACTAGTCCTTTGATGGATATGTGGTTTGTACATATTTTTTGTCAGTCTTTATCTTTTCATTTCATTCTCTTAATAGAGTCTTTCACAGAGCAAAAAGGTGTAAATTTTGATTAATCCAATAGATCAATGTTTTTCTTTTATAGGTTGTTCTTTTTTGTATTAAGTTTAACAAGAACACTTTGCCTAGCTTTAGATGCCAATGGTTTTTTTTGTGTGTGTGTGGAAAAATTTTATAGTTTTTATAATTTACATTTAAGTTCATGGTCCATTTTAAATTAATTTTCATATAAGCTAGTTGTTTTTTGTTGTTGTTTTCTCAATGAGTGTCCGTTTGCTTCGGGACCAAGAATGTCTTTTCTCTCTTAAATTGTTTTCTTTCACCTTTGTCAAAATTTAGTTAAGCATATTTGTGTGGGTCTATTTCTTTGTGGGTCTGTTTCTATTTCTGCTTCTAATTTTGTCCCATAGATCACTGTCTCTTCCTTCACCAGTACTGTACAATCTTGATTACTGTAGCTGTATAATAAGTCTTGAAATTGGGTAGACTGATTCTTCCCAGTTTAATCTTCTTTTAAAAAAATTTTTTTTTAGCTTTGTATTTCCTTCGCCTTTCCATAAAAATTTTTAGGATAATTTTGTCCGTATTTACTAAAAATCTTGCTGAGATTCTGATAGGAATTATGTTAAACCTCTCTATCAGTTTGAGAAGACTTGATGTCTTTACCATGTTAGTCTTTTAGTTCATGAACTTGGTATGTCTATCTGTGTTGTTAGATCTTTGATTTCTTTCATCAGTGTTGCGTAGTTTTCAGCATGCATTCTGTACATGTAGATTTACATCTGAGCATCTCAGTTTTTTTTTTTTTTTTTCAGTGATAGTAAATGAGCTCATTTATTAGCTCAAGGAGTTCTTTTTGTTTTCTGTTTTTTTGGGTAGTTTCCTTGAGATTTTACATAGACAGTCATATTTATCCTTTTTTTGAGATGGAGTCTCACTCTGTCACCCAGGCTGGAGTGCAGTGGCGCGATCTCGGCTCCCTGCAACCTCTGCTTCCCAGGTCCAAGCGATTCTCCTGCCTCAGCCTCCCAAGTAGCTGGGACTACAGGCATGCGCCACCATGCCCAGCTAATTTTTGCATTTTTTAGTAGAGACAAGGTTTCACCATGTTGGCCAGGCTGGTCTCGAACTCCTGACCTCAGGTGATCTACCCACCTTGGCCTCCTAAAATGCTGGAATTAGAGGCGTGAACTGCTGTGCTCAGCCTGTATTTTTATTTTCATTCCGTTGTGTGTATTTTAATTTCCCTTGAGACTCCCTCTTTTACCCATGGATTATTTTGAAGTATATTTTTGGGTTTCTAAGTATTTGGAAGATTTCTGTTATTTTTCTTTTATTGATTTCTAGCTTGATTCTTTGTGGTCAAAGAAACATACTATATGTGATTT	*	NM:i:2	AS:i:1095	XS:i:43	SA:Z:2,233376921,-,1091S13M1D994M,60,2;
NODE_1_length_2098_cov_6.69848_ID_3108	2064	2	233376921	60	1091H13M1D994M	*	0	0	CTATTTCTGCTTCTAATTTTGTCCCATAGATCACTGTCTCTTCCTTCACCAGTACTGTACAATCTTGATTACTGTAGCTGTATAATAAGTCTTGAAATTGGGTAGACTGATTCTTCCCAGTTTAATCTTCTTTTAAAAAAATTTTTTTTTAGCTTTGTATTTCCTTCGCCTTTCCATAAAAATTTTTAGGATAATTTTGTCCGTATTTACTAAAAATCTTGCTGAGATTCTGATAGGAATTATGTTAAACCTCTCTATCAGTTTGAGAAGACTTGATGTCTTTACCATGTTAGTCTTTTAGTTCATGAACTTGGTATGTCTATCTGTGTTGTTAGATCTTTGATTTCTTTCATCAGTGTTGCGTAGTTTTCAGCATGCATTCTGTACATGTAGATTTACATCTGAGCATCTCAGTTTTTTTTTTTTTTTTTCAGTGATAGTAAATGAGCTCATTTATTAGCTCAAGGAGTTCTTTTTGTTTTCTGTTTTTTTGGGTAGTTTCCTTGAGATTTTACATAGACAGTCATATTTATCCTTTTTTTGAGATGGAGTCTCACTCTGTCACCCAGGCTGGAGTGCAGTGGCGCGATCTCGGCTCCCTGCAACCTCTGCTTCCCAGGTCCAAGCGATTCTCCTGCCTCAGCCTCCCAAGTAGCTGGGACTACAGGCATGCGCCACCATGCCCAGCTAATTTTTGCATTTTTTAGTAGAGACAAGGTTTCACCATGTTGGCCAGGCTGGTCTCGAACTCCTGACCTCAGGTGATCTACCCACCTTGGCCTCCTAAAATGCTGGAATTAGAGGCGTGAACTGCTGTGCTCAGCCTGTATTTTTATTTTCATTCCGTTGTGTGTATTTTAATTTCCCTTGAGACTCCCTCTTTTACCCATGGATTATTTTGAAGTATATTTTTGGGTTTCTAAGTATTTGGAAGATTTCTGTTATTTTTCTTTTATTGATTTCTAGCTTGATTCTTTGTGGTCAAAGAAACATACTATATGTGATTT	*	NM:i:2	AS:i:995	XS:i:0	SA:Z:2,233375525,-,1105M993S,60,2;
"""
            sam_file_name = temp_dir + "/test.sam"
            sam_file = open(sam_file_name, "w")
            sam_file.write(sam_records)
            sam_file.close()
            print sam_file_name
            (samfile, records_by_name) = validate_local_assembly.read_sam_file(sam_file_name)
            
            deletion_intervals = validate_local_assembly.get_deletion_intervals(samfile, records_by_name)
            print deletion_intervals
            self.assertEqual(2, len(deletion_intervals))
            expected_deletion1 = pybedtools.Interval("2", 233376629, 233376920, "NODE_1_length_2155_cov_7.27857_ID_3018-1")
            expected_deletion2 = pybedtools.Interval("2", 233376933, 233376934, "NODE_1_length_2155_cov_7.27857_ID_3018-2")
            print deletion_intervals[0]
            print expected_deletion1
            self.assertEqual(expected_deletion1, deletion_intervals[0])
            print deletion_intervals[1]
            print expected_deletion2
            self.assertEqual(expected_deletion2, deletion_intervals[1])
        finally:
            shutil.rmtree(temp_dir)
    def test_deletion_between_alignments(self):
        # need to create temp files here until I can figure out how to mock out pysam
        temp_dir = tempfile.mkdtemp()
        
        try :
            sam_records = \
"""@SQ	SN:2	LN:242951149
NODE_1_length_2098_cov_6.69848_ID_3108	16	2	233375525	60	1105M993S	*	0	0	CATTCTAGGAATAAATAAATCCCACTTGGATTATGGTGTATAATCTTTTTAATATATTGCTGAATTTGCTTTGCTAGTATTTTGTTAAGAATTTTTGCATCAGTGTTCATAAAGGATATTGGTCTATATTTTTCTTTAGAGCCTTTGTCTGACTTTGATACCAGGGTAATTTTGGCTTCATAGAATTAGTTAGGAAGTGTAGGTAATCTTGCTTCCTTTTCAGTTTTTGGAAAAAGTTTGAGAAGGATTATGATTAGTTTTTCCTTAAATGTTTGGTGCGTTTTTCACCAGAGAAGCCATCAGGGCAGGGCTTTCCTTTGTCCAGAGATTCTTTTTTTTTTTTTTTAAAGACAGAGTCTAGCTCTGTCGGCCAGGCTGGAGTGCAGTGGCACAGTCTCAGCCCACTGCAACCGCTGTCTCCCGAGCTCAAGCAGTTCTCCTGCCTCAGCCTCCTGAGTAGCTGAGATTACAGGCGTGTGCCACCATGCCCGGCTCATTTTAGTAGAGATGGGGTTTCACCATGTTGGCCAGGCTGGTCTCGAACTCCTGACCTCAGATAATCTGCCTGCCTTGGCCTCCCAAACTGCTGGGATTACTGGCGTGAGCCACTGTGCCCGGCCTGTTCAGAGATTTTTGATTACTGATTCTAGATACTAGTCCTTTGATGGATATGTGGTTTGTACATATTTTTTGTCAGTCTTTATCTTTTCATTTCATTCTCTTAATAGAGTCTTTCACAGAGCAAAAAGGTGTAAATTTTGATTAATCCAATAGATCAATGTTTTTCTTTTATAGGTTGTTCTTTTTTGTATTAAGTTTAACAAGAACACTTTGCCTAGCTTTAGATGCCAATGGTTTTTTTTGTGTGTGTGTGGAAAAATTTTATAGTTTTTATAATTTACATTTAAGTTCATGGTCCATTTTAAATTAATTTTCATATAAGCTAGTTGTTTTTTGTTGTTGTTTTCTCAATGAGTGTCCGTTTGCTTCGGGACCAAGAATGTCTTTTCTCTCTTAAATTGTTTTCTTTCACCTTTGTCAAAATTTAGTTAAGCATATTTGTGTGGGTCTATTTCTTTGTGGGTCTGTTTCTATTTCTGCTTCTAATTTTGTCCCATAGATCACTGTCTCTTCCTTCACCAGTACTGTACAATCTTGATTACTGTAGCTGTATAATAAGTCTTGAAATTGGGTAGACTGATTCTTCCCAGTTTAATCTTCTTTTAAAAAAATTTTTTTTTAGCTTTGTATTTCCTTCGCCTTTCCATAAAAATTTTTAGGATAATTTTGTCCGTATTTACTAAAAATCTTGCTGAGATTCTGATAGGAATTATGTTAAACCTCTCTATCAGTTTGAGAAGACTTGATGTCTTTACCATGTTAGTCTTTTAGTTCATGAACTTGGTATGTCTATCTGTGTTGTTAGATCTTTGATTTCTTTCATCAGTGTTGCGTAGTTTTCAGCATGCATTCTGTACATGTAGATTTACATCTGAGCATCTCAGTTTTTTTTTTTTTTTTTCAGTGATAGTAAATGAGCTCATTTATTAGCTCAAGGAGTTCTTTTTGTTTTCTGTTTTTTTGGGTAGTTTCCTTGAGATTTTACATAGACAGTCATATTTATCCTTTTTTTGAGATGGAGTCTCACTCTGTCACCCAGGCTGGAGTGCAGTGGCGCGATCTCGGCTCCCTGCAACCTCTGCTTCCCAGGTCCAAGCGATTCTCCTGCCTCAGCCTCCCAAGTAGCTGGGACTACAGGCATGCGCCACCATGCCCAGCTAATTTTTGCATTTTTTAGTAGAGACAAGGTTTCACCATGTTGGCCAGGCTGGTCTCGAACTCCTGACCTCAGGTGATCTACCCACCTTGGCCTCCTAAAATGCTGGAATTAGAGGCGTGAACTGCTGTGCTCAGCCTGTATTTTTATTTTCATTCCGTTGTGTGTATTTTAATTTCCCTTGAGACTCCCTCTTTTACCCATGGATTATTTTGAAGTATATTTTTGGGTTTCTAAGTATTTGGAAGATTTCTGTTATTTTTCTTTTATTGATTTCTAGCTTGATTCTTTGTGGTCAAAGAAACATACTATATGTGATTT	*	NM:i:2	AS:i:1095	XS:i:43	SA:Z:2,233376921,-,1091S13M1D994M,60,2;
NODE_1_length_2098_cov_6.69848_ID_3108	2064	2	233376921	60	1091H13M1D994M	*	0	0	CTATTTCTGCTTCTAATTTTGTCCCATAGATCACTGTCTCTTCCTTCACCAGTACTGTACAATCTTGATTACTGTAGCTGTATAATAAGTCTTGAAATTGGGTAGACTGATTCTTCCCAGTTTAATCTTCTTTTAAAAAAATTTTTTTTTAGCTTTGTATTTCCTTCGCCTTTCCATAAAAATTTTTAGGATAATTTTGTCCGTATTTACTAAAAATCTTGCTGAGATTCTGATAGGAATTATGTTAAACCTCTCTATCAGTTTGAGAAGACTTGATGTCTTTACCATGTTAGTCTTTTAGTTCATGAACTTGGTATGTCTATCTGTGTTGTTAGATCTTTGATTTCTTTCATCAGTGTTGCGTAGTTTTCAGCATGCATTCTGTACATGTAGATTTACATCTGAGCATCTCAGTTTTTTTTTTTTTTTTTCAGTGATAGTAAATGAGCTCATTTATTAGCTCAAGGAGTTCTTTTTGTTTTCTGTTTTTTTGGGTAGTTTCCTTGAGATTTTACATAGACAGTCATATTTATCCTTTTTTTGAGATGGAGTCTCACTCTGTCACCCAGGCTGGAGTGCAGTGGCGCGATCTCGGCTCCCTGCAACCTCTGCTTCCCAGGTCCAAGCGATTCTCCTGCCTCAGCCTCCCAAGTAGCTGGGACTACAGGCATGCGCCACCATGCCCAGCTAATTTTTGCATTTTTTAGTAGAGACAAGGTTTCACCATGTTGGCCAGGCTGGTCTCGAACTCCTGACCTCAGGTGATCTACCCACCTTGGCCTCCTAAAATGCTGGAATTAGAGGCGTGAACTGCTGTGCTCAGCCTGTATTTTTATTTTCATTCCGTTGTGTGTATTTTAATTTCCCTTGAGACTCCCTCTTTTACCCATGGATTATTTTGAAGTATATTTTTGGGTTTCTAAGTATTTGGAAGATTTCTGTTATTTTTCTTTTATTGATTTCTAGCTTGATTCTTTGTGGTCAAAGAAACATACTATATGTGATTT	*	NM:i:2	AS:i:995	XS:i:0	SA:Z:2,233375525,-,1105M993S,60,2;
"""
            sam_file_name = temp_dir + "/test.sam"
            sam_file = open(sam_file_name, "w")
            sam_file.write(sam_records)
            sam_file.close()
            print sam_file_name
            (samfile, records_by_name) = validate_local_assembly.read_sam_file(sam_file_name)
            
            deletion_intervals = validate_local_assembly.get_deletion_intervals(samfile, records_by_name)
            print deletion_intervals
            self.assertEqual(2, len(deletion_intervals))
            expected_deletion1 = pybedtools.Interval("2", 233376629, 233376920, "NODE_1_length_2155_cov_7.27857_ID_3018-1")
            expected_deletion2 = pybedtools.Interval("2", 233376933, 233376934, "NODE_1_length_2155_cov_7.27857_ID_3018-2")
            print deletion_intervals[0]
            print expected_deletion1
            self.assertEqual(expected_deletion1, deletion_intervals[0])
            print deletion_intervals[1]
            print expected_deletion2
            self.assertEqual(expected_deletion2, deletion_intervals[1])
        finally:
            shutil.rmtree(temp_dir)

    def test_deletion_between_alignments(self):
        # need to create temp files here until I can figure out how to mock out pysam
        temp_dir = tempfile.mkdtemp()
        
        try :
            sam_records = \
"""@SQ	SN:2	LN:242951149
NODE_1_length_2098_cov_6.69848_ID_3108	16	2	233375525	60	1105M993S	*	0	0	CATTCTAGGAATAAATAAATCCCACTTGGATTATGGTGTATAATCTTTTTAATATATTGCTGAATTTGCTTTGCTAGTATTTTGTTAAGAATTTTTGCATCAGTGTTCATAAAGGATATTGGTCTATATTTTTCTTTAGAGCCTTTGTCTGACTTTGATACCAGGGTAATTTTGGCTTCATAGAATTAGTTAGGAAGTGTAGGTAATCTTGCTTCCTTTTCAGTTTTTGGAAAAAGTTTGAGAAGGATTATGATTAGTTTTTCCTTAAATGTTTGGTGCGTTTTTCACCAGAGAAGCCATCAGGGCAGGGCTTTCCTTTGTCCAGAGATTCTTTTTTTTTTTTTTTAAAGACAGAGTCTAGCTCTGTCGGCCAGGCTGGAGTGCAGTGGCACAGTCTCAGCCCACTGCAACCGCTGTCTCCCGAGCTCAAGCAGTTCTCCTGCCTCAGCCTCCTGAGTAGCTGAGATTACAGGCGTGTGCCACCATGCCCGGCTCATTTTAGTAGAGATGGGGTTTCACCATGTTGGCCAGGCTGGTCTCGAACTCCTGACCTCAGATAATCTGCCTGCCTTGGCCTCCCAAACTGCTGGGATTACTGGCGTGAGCCACTGTGCCCGGCCTGTTCAGAGATTTTTGATTACTGATTCTAGATACTAGTCCTTTGATGGATATGTGGTTTGTACATATTTTTTGTCAGTCTTTATCTTTTCATTTCATTCTCTTAATAGAGTCTTTCACAGAGCAAAAAGGTGTAAATTTTGATTAATCCAATAGATCAATGTTTTTCTTTTATAGGTTGTTCTTTTTTGTATTAAGTTTAACAAGAACACTTTGCCTAGCTTTAGATGCCAATGGTTTTTTTTGTGTGTGTGTGGAAAAATTTTATAGTTTTTATAATTTACATTTAAGTTCATGGTCCATTTTAAATTAATTTTCATATAAGCTAGTTGTTTTTTGTTGTTGTTTTCTCAATGAGTGTCCGTTTGCTTCGGGACCAAGAATGTCTTTTCTCTCTTAAATTGTTTTCTTTCACCTTTGTCAAAATTTAGTTAAGCATATTTGTGTGGGTCTATTTCTTTGTGGGTCTGTTTCTATTTCTGCTTCTAATTTTGTCCCATAGATCACTGTCTCTTCCTTCACCAGTACTGTACAATCTTGATTACTGTAGCTGTATAATAAGTCTTGAAATTGGGTAGACTGATTCTTCCCAGTTTAATCTTCTTTTAAAAAAATTTTTTTTTAGCTTTGTATTTCCTTCGCCTTTCCATAAAAATTTTTAGGATAATTTTGTCCGTATTTACTAAAAATCTTGCTGAGATTCTGATAGGAATTATGTTAAACCTCTCTATCAGTTTGAGAAGACTTGATGTCTTTACCATGTTAGTCTTTTAGTTCATGAACTTGGTATGTCTATCTGTGTTGTTAGATCTTTGATTTCTTTCATCAGTGTTGCGTAGTTTTCAGCATGCATTCTGTACATGTAGATTTACATCTGAGCATCTCAGTTTTTTTTTTTTTTTTTCAGTGATAGTAAATGAGCTCATTTATTAGCTCAAGGAGTTCTTTTTGTTTTCTGTTTTTTTGGGTAGTTTCCTTGAGATTTTACATAGACAGTCATATTTATCCTTTTTTTGAGATGGAGTCTCACTCTGTCACCCAGGCTGGAGTGCAGTGGCGCGATCTCGGCTCCCTGCAACCTCTGCTTCCCAGGTCCAAGCGATTCTCCTGCCTCAGCCTCCCAAGTAGCTGGGACTACAGGCATGCGCCACCATGCCCAGCTAATTTTTGCATTTTTTAGTAGAGACAAGGTTTCACCATGTTGGCCAGGCTGGTCTCGAACTCCTGACCTCAGGTGATCTACCCACCTTGGCCTCCTAAAATGCTGGAATTAGAGGCGTGAACTGCTGTGCTCAGCCTGTATTTTTATTTTCATTCCGTTGTGTGTATTTTAATTTCCCTTGAGACTCCCTCTTTTACCCATGGATTATTTTGAAGTATATTTTTGGGTTTCTAAGTATTTGGAAGATTTCTGTTATTTTTCTTTTATTGATTTCTAGCTTGATTCTTTGTGGTCAAAGAAACATACTATATGTGATTT	*	NM:i:2	AS:i:1095	XS:i:43	SA:Z:2,233376921,-,1091S13M1D994M,60,2;
NODE_1_length_2098_cov_6.69848_ID_3108	2064	2	233376921	60	1091H13M1D994M	*	0	0	CTATTTCTGCTTCTAATTTTGTCCCATAGATCACTGTCTCTTCCTTCACCAGTACTGTACAATCTTGATTACTGTAGCTGTATAATAAGTCTTGAAATTGGGTAGACTGATTCTTCCCAGTTTAATCTTCTTTTAAAAAAATTTTTTTTTAGCTTTGTATTTCCTTCGCCTTTCCATAAAAATTTTTAGGATAATTTTGTCCGTATTTACTAAAAATCTTGCTGAGATTCTGATAGGAATTATGTTAAACCTCTCTATCAGTTTGAGAAGACTTGATGTCTTTACCATGTTAGTCTTTTAGTTCATGAACTTGGTATGTCTATCTGTGTTGTTAGATCTTTGATTTCTTTCATCAGTGTTGCGTAGTTTTCAGCATGCATTCTGTACATGTAGATTTACATCTGAGCATCTCAGTTTTTTTTTTTTTTTTTCAGTGATAGTAAATGAGCTCATTTATTAGCTCAAGGAGTTCTTTTTGTTTTCTGTTTTTTTGGGTAGTTTCCTTGAGATTTTACATAGACAGTCATATTTATCCTTTTTTTGAGATGGAGTCTCACTCTGTCACCCAGGCTGGAGTGCAGTGGCGCGATCTCGGCTCCCTGCAACCTCTGCTTCCCAGGTCCAAGCGATTCTCCTGCCTCAGCCTCCCAAGTAGCTGGGACTACAGGCATGCGCCACCATGCCCAGCTAATTTTTGCATTTTTTAGTAGAGACAAGGTTTCACCATGTTGGCCAGGCTGGTCTCGAACTCCTGACCTCAGGTGATCTACCCACCTTGGCCTCCTAAAATGCTGGAATTAGAGGCGTGAACTGCTGTGCTCAGCCTGTATTTTTATTTTCATTCCGTTGTGTGTATTTTAATTTCCCTTGAGACTCCCTCTTTTACCCATGGATTATTTTGAAGTATATTTTTGGGTTTCTAAGTATTTGGAAGATTTCTGTTATTTTTCTTTTATTGATTTCTAGCTTGATTCTTTGTGGTCAAAGAAACATACTATATGTGATTT	*	NM:i:2	AS:i:995	XS:i:0	SA:Z:2,233375525,-,1105M993S,60,2;
"""
            sam_file_name = temp_dir + "/test.sam"
            sam_file = open(sam_file_name, "w")
            sam_file.write(sam_records)
            sam_file.close()
            print sam_file_name
            (samfile, records_by_name) = validate_local_assembly.read_sam_file(sam_file_name)
            
            deletion_intervals = validate_local_assembly.get_deletion_intervals(samfile, records_by_name)
            print deletion_intervals
            self.assertEqual(2, len(deletion_intervals))
            expected_deletion1 = pybedtools.Interval("2", 233376629, 233376920, "NODE_1_length_2155_cov_7.27857_ID_3018-1")
            expected_deletion2 = pybedtools.Interval("2", 233376933, 233376934, "NODE_1_length_2155_cov_7.27857_ID_3018-2")
            print deletion_intervals[0]
            print expected_deletion1
            self.assertEqual(expected_deletion1, deletion_intervals[0])
            print deletion_intervals[1]
            print expected_deletion2
            self.assertEqual(expected_deletion2, deletion_intervals[1])
        finally:
            shutil.rmtree(temp_dir)

    def test_weird_palindrome(self):
        # need to create temp files here until I can figure out how to mock out pysam
        temp_dir = tempfile.mkdtemp()
        
        try :
            sam_records = \
"""@SQ	SN:2	LN:242951149
NODE_1_length_1777_cov_10.6794_ID_3450	0	2	83575057	60	1148M629S	*	0	0	CCAGGAATGAATACTTGGAAATACAAATTGAAAACACAATATTATTTATATTAACACCAAAAAAATTATATCATTAGGTTTAAGCGTAGCAAAGTATATAAAAGATCTATATGAGATAAGCTACAGAAACCCAATGATATAAATCAAATAAGATCTAAATAAACAGATATCCCATGGTTATGAATAGATATTGTTAAGATAACCTTTATTTTCCACTTGATCTATAGATTCAATATAATACCAATCAAAATCCAAGCAAATTATTTTGTGGATATTGAAAGCTAATGCTAATAATTATGAAGATTAGTAAAAAGACTCAGAATAGCCAACTGAATTCTGAAGAAGAACAAAGTTGAAAAACTGACACTACCCAACTTTAACCTATTACTATGGATTACAAGTAAATCTATAGTAATCAAGACAGTCTGGTATTGGTGAAAGAATATACAAATAATTCAGTGGAACCAAATTGATGCCCAGGAATAGACCAACAAAAATACAGTCAATTGATCTTTGACACAAAGCAAAGGCAATTAAATGAAGAAAGAAGTCTTTTCAGTAAATGGTGCTGCAACAATCACATATCCAAGTACCAAAAAAATAAATTTTGGCACAGAACTTATGTCCTAAAAATTAAATCATGGTGTATATGTGCCACATTTTCTTAATCCAGTCTATCATTGTTGGACATTTGGGTTGGTTCCAAGTCTTTGCTAATTGTGAGTAGTGCTGCAATAAAAATACATGTGCATGTGTCTTTATAGCAGCATGATTTGTATTCCTTTTGGCACAATTACACCATGGAATACTATGCAACCATAAAACATAATGAGTTCATGCCCTTTGTAGGGACATGGATGAAGCTGGAAACCATCATTCTCAGCAAACTATAGCAAGGACAAAAAACCAAACACCGCATGTTCTCACTCATAGGTGGGAATTGAACAATGAGAACACTTGGACACAGGAAGGGGAACATCACATGCCGGAACCTGTTGTAGGGTGGGAGGAGGAGGGAGGGATAGCATTAGGAGATATACCTAATGTAAATGATGAGTTAATGGGTGCAGCACACCAACATAGCACATGTATACATATGTAACAAACCTGCCTGTTGTGCACATGTACCGTAGAACTTAAAGTATAATCTATATATATATAATATATATATATAAAAGAAGCATCTGTATTATCTTGCAAAAGTCTTAAAATTTATTTTACATTTAGCTTTAAGATCTATTTTGATTTAATTTTTTAATATTTTACTTTTTAAGAAAAAACACATAATAATTGCATAAAAACATAAGTATGCAGTTTAATGAATTACTATAATAGGAACACTGTATAATCACCAATCAAGGTCAAGAAGTAGAATATTGGTAGGGCTCAGAAGCACCCCTACCTCCATGCCCATCCCCATCCACAGCTTCCTTTCTCTTCTTAAGCGTAGTCCCTTCCTGACATTTATGGGTAATGTGTAGTCACTTTTTTGTTGTTTATTCTCTATTTTTTTATTTTTATTTTTATTTTTTGGAGTTAGGGGACCTAAATTCTATATATATATAATATATATATAAAATATATATATAATATATATATAAAATATATATATAATATATATATAAAATATATATATAATATATATATAAAATATATATATAATATATATATAAAATATATATATAATATATATATAAAATATATATATAATATATATATAAAATATATATATAATATATATATAAAATATATATATAATATATATATAAAATATATA	*	NM:i:7	AS:i:1113	XS:i:186	SA:Z:2,83576358,-,718M1D265M794S,11,30;
NODE_1_length_1777_cov_10.6794_ID_3450	2064	2	83576358	11	718M1D265M794H	*	0	0	TATATATTTTATATATATATTATATATATATTTTATATATATATTATATATATATTTTATATATATATTATATATATATTTTATATATATATTATATATATATTTTATATATATATTATATATATATTTTATATATATATTATATATATATTTTATATATATATTATATATATATTTTATATATATATTATATATATATTTTATATATATATTATATATATATAGAATTTAGGTCCCCTAACTCCAAAAAATAAAAATAAAAATAAAAAAATAGAGAATAAACAACAAAAAAGTGACTACACATTACCCATAAATGTCAGGAAGGGACTACGCTTAAGAAGAGAAAGGAAGCTGTGGATGGGGATGGGCATGGAGGTAGGGGTGCTTCTGAGCCCTACCAATATTCTACTTCTTGACCTTGATTGGTGATTATACAGTGTTCCTATTATAGTAATTCATTAAACTGCATACTTATGTTTTTATGCAATTATTATGTGTTTTTTCTTAAAAAGTAAAATATTAAAAAATTAAATCAAAATAGATCTTAAAGCTAAATGTAAAATAAATTTTAAGACTTTTGCAAGATAATACAGATGCTTCTTTTATATATATATATTATATATATATAGATTATACTTTAAGTTCTACGGTACATGTGCACAACAGGCAGGTTTGTTACATATGTATACATGTGCTATGTTGGTGTGCTGCACCCATTAACTCATCATTTACATTAGGTATATCTCCTAATGCTATCCCTCCCTCCTCCTCCCACCCTACAACAGGTTCCGGCATGTGATGTTCCCCTTCCTGTGTCCAAGTGTTCTCATTGTTCAATTCCCACCTATGAGTGAGAACATGCGGTGTTTGGTTTTTTGTCCTTGCTATAGTTTGCTGAGAATGATGGTTTCCAGCTTCATCCATGTCCCTACAAAGGGCATGAACTCATTATGTTTTATGGTTGCATAGTATTCCATGGTGTA	*	NM:i:30	AS:i:831	XS:i:787	SA:Z:2,83575057,+,1148M629S,60,7;
"""
            sam_file_name = temp_dir + "/test.sam"
            sam_file = open(sam_file_name, "w")
            sam_file.write(sam_records)
            sam_file.close()
            print sam_file_name
            (samfile, records_by_name) = validate_local_assembly.read_sam_file(sam_file_name)
            
            deletion_intervals = validate_local_assembly.get_deletion_intervals(samfile, records_by_name)
            print deletion_intervals
            self.assertEqual(1, len(deletion_intervals))
            expected_deletion1 = pybedtools.Interval("2", 83577075, 83577076, "NODE_1_length_2155_cov_7.27857_ID_3018-1")
            print deletion_intervals[0]
            print expected_deletion1
            self.assertEqual(expected_deletion1, deletion_intervals[0])
        finally:
            shutil.rmtree(temp_dir)

    def test_overlapping_alignments(self):
        # need to create temp files here until I can figure out how to mock out pysam
        temp_dir = tempfile.mkdtemp()
        
        try :
            sam_records = \
"""@SQ	SN:2	LN:242951149
NODE_1_length_2059_cov_6.98902_ID_3122	16	2	1042120	60	625S16M1D35M1D35M1D16M4D32M4D62M5I38M4I176M1D80M4I931M	*	0	0	AGTAGCAGGGACAATAATTTATCAGCTTTATTATGAAAATGTTAGTCATTCATAGAGATCAAACCATGTTGCAGAAATTCAGATTTGGAAGAAACCCGGAACAGTAATTGATCTTATCCATTCTCCTAGGGTTATCATGAAGGATTTAACAGTGTGAACATTTCCATTTTCTAAGGTGACCGTGCCTTTATTTTGATGGCTCTGAATTGTATTTTCCATTCAGTTGCAGCACCAGCCCTCCCTTGATGTGCAGAATGAGCAGGTGTCTTTAGGAGTGAGGCTGGGGCACGACTGAATGCAAATTAGCCCTTGCCCAGACAGAAGACAGTCTCATCTCACCAGCGGTATTTCCTGGAGGGTCTGATCTGTCCTGGACATGCGGTGGGTCGCGGGGACTGCAGCTTCATCTCTCAGACGGCTCCTTGTCCCCACTGTTTCAGAGTGGACACAGGGAGAAGGGGTGAGGCAGCTGGGACCCAGGCGGCTCCTCTGGCTCCCATCTTGGTGGGAAGCTGGTCCACACCACCAAGAGCACGGGGAGGACTCCCTGCGGATGGGCTCCAGCTTCATTCTGAAGGAAGGAGTTGAATAAGTTAAGGAAAGTGGATGAGAGGGAGAGTGCATGGGGGAGGGAGGGAGAGGCCGCACTGTGCTGTGGGGAGGGAGGGAGGGAGAGGCCGCGCCGTGCTGTGGGGAGGGAGGGAGGGAGAGGCCGCGCCGTGCTGTGGGGAGGGAGGGAGAGGGCGGCGCCCCGCTGTGGGGAGGGAGGGAGAGGGCGGCGCCCCGCTGTGGGGAGGGAGGGAGGGAGAGCGCAGCGCCACGCTGTGCTGTGGGGAGGGAGGGAGAGCGTCGCGCTGTGCTGTGGGGAGGGAGGGAGGGAGAGCGCAGCGCCACGCTGTGGGGAGGGAGGGAGGGAGAGCACGGCGCCCCGCTGTGGGGAGGGAGGGAGGGAGAGCGCAGCGCCACGCTGTGGGGAGGGAGGGAGAGCGTCGCGCTGTGCTGTGGGGCGGGAGGGAGAGGGCGGCGCCACGCTGTGGGGAGGGAAGAGGGAGAGCACGGCGCCCCGCTGTGCTGTGGGGAGGGAGGGAGAGCGTCGCGCTGTGCTGTGGGGAGGGAGGGAGGGAGGAGGAGCGCGGCGCCCCGCTGTGGGGAGGGAGGGAGGGAGAGCGCGGGGCCACCCCGTGGGGAGGGAGGAGAGGGCGGCGCAGCGCTGTGCTGCGGGGAACGATGCATGCCGCTTCCTCCAGGACCTGTGAGTGCAGGGTTGGACGGGGCTTTCTGAGCTTGGGAGCTGGGGCCAGGCCTGTCCTGCAAGCAGCTGGCCTATTTCTGCTAGAAGCTGTGAATTACAACTGTAATGTATTTATTCTAAGGACTAACATGGCACTATTCATAAAAGACTGTGTTAATATTTGCACTGTGCTTTCTGTAATTATACCTCATCAACTTCATTTTTATCTCTATTCTGGCTGTATTTCATATTCCTGATTCTCTTCTGTATGCCTGAAAGGAGAGTTGGTGTATTAGTCAGGGTTCTCCAGAGGGACAGAATTAATAGGATAGACACATGCATGTAAAGGGGAGCATATTAAGTAGTAGTAACTCACACGTTCACAAGTGTGAACATAGTCCGTCTGCAAGCTGAGGAGCAAGGAAGCCAGTTCGAGTCCCAAAGCTGGAGTCCGATGTGCGAGGGCAGGAAGCATCCAGTGCGGGAGAAAGATGCTAAGCCATTCTAGCCTCTTTATGTTTTTCTGTCTACTTTATATTTTGGCTGCACTGGCCCTGATTAGATGGTGCCCACTGAAATTAAGGGTGGGTCTGTTTTTCCCAGTCCACTGACTCAAATGTTAGTCTCCTTTGGCAACACCCTCACAGACACACCCAGGATCAGTACTTTGCATCCTTCAGTCCAATCAAGTTGACACTCACTATTAACCATCACAAGTGGCAAAATATAAATGCCTCTAAATGTGCACAGGATCTTTCTATATTTAAACAATTAGTATAAAGTTGCATTTCAAAAAAATCTTCCAGCTGGACAACATGGTGAGACCCC	*	NM:i:43	AS:i:1252	XS:i:133	SA:Z:2,1041289,-,972M4D32M4D38M3I8M5I66M4I66M865S,4,49;
NODE_1_length_2059_cov_6.98902_ID_3122	2064	2	1041289	4	972M4D32M4D38M3I8M5I66M4I66M865H	*	0	0	AGTAGCAGGGACAATAATTTATCAGCTTTATTATGAAAATGTTAGTCATTCATAGAGATCAAACCATGTTGCAGAAATTCAGATTTGGAAGAAACCCGGAACAGTAATTGATCTTATCCATTCTCCTAGGGTTATCATGAAGGATTTAACAGTGTGAACATTTCCATTTTCTAAGGTGACCGTGCCTTTATTTTGATGGCTCTGAATTGTATTTTCCATTCAGTTGCAGCACCAGCCCTCCCTTGATGTGCAGAATGAGCAGGTGTCTTTAGGAGTGAGGCTGGGGCACGACTGAATGCAAATTAGCCCTTGCCCAGACAGAAGACAGTCTCATCTCACCAGCGGTATTTCCTGGAGGGTCTGATCTGTCCTGGACATGCGGTGGGTCGCGGGGACTGCAGCTTCATCTCTCAGACGGCTCCTTGTCCCCACTGTTTCAGAGTGGACACAGGGAGAAGGGGTGAGGCAGCTGGGACCCAGGCGGCTCCTCTGGCTCCCATCTTGGTGGGAAGCTGGTCCACACCACCAAGAGCACGGGGAGGACTCCCTGCGGATGGGCTCCAGCTTCATTCTGAAGGAAGGAGTTGAATAAGTTAAGGAAAGTGGATGAGAGGGAGAGTGCATGGGGGAGGGAGGGAGAGGCCGCACTGTGCTGTGGGGAGGGAGGGAGGGAGAGGCCGCGCCGTGCTGTGGGGAGGGAGGGAGGGAGAGGCCGCGCCGTGCTGTGGGGAGGGAGGGAGAGGGCGGCGCCCCGCTGTGGGGAGGGAGGGAGAGGGCGGCGCCCCGCTGTGGGGAGGGAGGGAGGGAGAGCGCAGCGCCACGCTGTGCTGTGGGGAGGGAGGGAGAGCGTCGCGCTGTGCTGTGGGGAGGGAGGGAGGGAGAGCGCAGCGCCACGCTGTGGGGAGGGAGGGAGGGAGAGCACGGCGCCCCGCTGTGGGGAGGGAGGGAGGGAGAGCGCAGCGCCACGCTGTGGGGAGGGAGGGAGAGCGTCGCGCTGTGCTGTGGGGCGGGAGGGAGAGGGCGGCGCCACGCTGTGGGGAGGGAAGAGGGAGAGCACGGCGCCCCGCTGTGCTGTGGGGAGGGAGGGAGAGCGTCGCGCTGTGCTGTGGGGAGGGAGGGAGGGAGGAGGAGCGCGGCGCCCCGCTGTGGGGAGGGAGGGAGGGAGAGCGCGGGGCCACCCCGTGGGGAGGGAGG	*	NM:i:49	AS:i:987	XS:i:969	SA:Z:2,1042120,-,625S16M1D35M1D35M1D16M4D32M4D62M5I38M4I176M1D80M4I931M,60,43;
"""
            sam_file_name = temp_dir + "/test.sam"
            sam_file = open(sam_file_name, "w")
            sam_file.write(sam_records)
            sam_file.close()
            print sam_file_name
            (samfile, records_by_name) = validate_local_assembly.read_sam_file(sam_file_name)
            
            deletion_intervals = validate_local_assembly.get_deletion_intervals(samfile, records_by_name)
            print deletion_intervals
            self.assertEqual(2, len(deletion_intervals))
            expected_deletion1 = pybedtools.Interval("2", 1042260, 1042264, "NODE_1_length_2059_cov_6.98902_ID_3122-1")
            print deletion_intervals[0]
            print expected_deletion1
            self.assertEqual(expected_deletion1, deletion_intervals[0])
            expected_deletion2 = pybedtools.Interval("2", 1042296, 1042300, "NODE_1_length_2059_cov_6.98902_ID_3122-2")
            self.assertEqual(expected_deletion2, deletion_intervals[1])
        finally:
            shutil.rmtree(temp_dir)

if __name__ == '__main__':
    unittest.main()
        
