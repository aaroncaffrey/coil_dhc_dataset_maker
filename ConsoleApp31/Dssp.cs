using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace Coil_DHC_Dataset_Maker
{
    public static class Dssp
    {
        public class DsspRecord
        {
            public const char dssp_default_secondary_structure = 'C';

            public string SequentualResidueNumber;
            public string PdbResidueSequenceIndex;
            public char iCode = ' ';
            public char Chain;
            public char AminoAcid;
            public char SecondaryStructure = dssp_default_secondary_structure;
            public char TurnsHelix3;
            public char TurnsHelix4;
            public char TurnsHelix5;
            public char GeometricalBend;
            public char Chirality;
            public string BetaBridgeLabel2;
            public string BridgePartner1;
            public string BridgePartner2;
            public char BetaSheetLabel;
            public string Acc;
            public string N_H___O1;
            public string O___H_N1;
            public string N_H___O2;
            public string O___H_N2;
            public string Tco;
            public string Kappa;
            public string Alpha;
            public string PHI;
            public string PSI;
            public string X_CA;
            public string Y_CA;
            public string Z_CA;

            public DsspRecord(string line)
            {
                SequentualResidueNumber = (1 > 0 && 5 > 0 && line.Length >= 1) ? line.Substring(1 - 1, (line.Length >= 5 ? ((5 - 1) + 1) : line.Length - (1 - 1))).Trim() : "";
                PdbResidueSequenceIndex = (6 > 0 && 10 > 0 && line.Length >= 6) ? line.Substring(6 - 1, (line.Length >= 10 ? ((10 - 6) + 1) : line.Length - (6 - 1))).Trim() : "";
                iCode = (11 > 0 && 11 > 0 && line.Length >= 11) ? line.Substring(11 - 1, (line.Length >= 11 ? ((11 - 11) + 1) : line.Length - (11 - 1)))[0] : ' ';
                Chain = (12 > 0 && 12 > 0 && line.Length >= 12) ? line.Substring(12 - 1, (line.Length >= 12 ? ((12 - 12) + 1) : line.Length - (12 - 1)))[0] : ' ';
                AminoAcid = (14 > 0 && 14 > 0 && line.Length >= 14) ? line.Substring(14 - 1, (line.Length >= 14 ? ((14 - 14) + 1) : line.Length - (14 - 1)))[0] : ' ';
                SecondaryStructure = (17 > 0 && 17 > 0 && line.Length >= 17) ? line.Substring(17 - 1, (line.Length >= 17 ? ((17 - 17) + 1) : line.Length - (17 - 1)))[0] : dssp_default_secondary_structure;
                TurnsHelix3 = (19 > 0 && 19 > 0 && line.Length >= 19) ? line.Substring(19 - 1, (line.Length >= 19 ? ((19 - 19) + 1) : line.Length - (19 - 1)))[0] : ' ';
                TurnsHelix4 = (20 > 0 && 20 > 0 && line.Length >= 20) ? line.Substring(20 - 1, (line.Length >= 20 ? ((20 - 20) + 1) : line.Length - (20 - 1)))[0] : ' ';
                TurnsHelix5 = (21 > 0 && 21 > 0 && line.Length >= 21) ? line.Substring(21 - 1, (line.Length >= 21 ? ((21 - 21) + 1) : line.Length - (21 - 1)))[0] : ' ';
                GeometricalBend = (22 > 0 && 22 > 0 && line.Length >= 22) ? line.Substring(22 - 1, (line.Length >= 22 ? ((22 - 22) + 1) : line.Length - (22 - 1)))[0] : ' ';
                Chirality = (23 > 0 && 23 > 0 && line.Length >= 23) ? line.Substring(23 - 1, (line.Length >= 23 ? ((23 - 23) + 1) : line.Length - (23 - 1)))[0] : ' ';
                BetaBridgeLabel2 = (24 > 0 && 25 > 0 && line.Length >= 24) ? line.Substring(24 - 1, (line.Length >= 25 ? ((25 - 24) + 1) : line.Length - (24 - 1))).Trim() : "";
                BridgePartner1 = (27 > 0 && 29 > 0 && line.Length >= 27) ? line.Substring(27 - 1, (line.Length >= 29 ? ((29 - 27) + 1) : line.Length - (27 - 1))).Trim() : "";
                BridgePartner2 = (31 > 0 && 33 > 0 && line.Length >= 31) ? line.Substring(31 - 1, (line.Length >= 33 ? ((33 - 31) + 1) : line.Length - (31 - 1))).Trim() : "";
                BetaSheetLabel = (34 > 0 && 34 > 0 && line.Length >= 34) ? line.Substring(34 - 1, (line.Length >= 34 ? ((34 - 34) + 1) : line.Length - (34 - 1)))[0] : ' ';
                Acc = (36 > 0 && 38 > 0 && line.Length >= 36) ? line.Substring(36 - 1, (line.Length >= 38 ? ((38 - 36) + 1) : line.Length - (36 - 1))).Trim() : "";
                N_H___O1 = (40 > 0 && 50 > 0 && line.Length >= 40) ? line.Substring(40 - 1, (line.Length >= 50 ? ((50 - 40) + 1) : line.Length - (40 - 1))).Trim() : "";
                O___H_N1 = (52 > 0 && 62 > 0 && line.Length >= 52) ? line.Substring(52 - 1, (line.Length >= 62 ? ((62 - 52) + 1) : line.Length - (52 - 1))).Trim() : "";
                N_H___O2 = (64 > 0 && 74 > 0 && line.Length >= 64) ? line.Substring(64 - 1, (line.Length >= 74 ? ((74 - 64) + 1) : line.Length - (64 - 1))).Trim() : "";
                O___H_N2 = (76 > 0 && 86 > 0 && line.Length >= 76) ? line.Substring(76 - 1, (line.Length >= 86 ? ((86 - 76) + 1) : line.Length - (76 - 1))).Trim() : "";
                Tco = (86 > 0 && 91 > 0 && line.Length >= 86) ? line.Substring(86 - 1, (line.Length >= 91 ? ((91 - 86) + 1) : line.Length - (86 - 1))).Trim() : "";
                Kappa = (93 > 0 && 97 > 0 && line.Length >= 93) ? line.Substring(93 - 1, (line.Length >= 97 ? ((97 - 93) + 1) : line.Length - (93 - 1))).Trim() : "";
                Alpha = (99 > 0 && 103 > 0 && line.Length >= 99) ? line.Substring(99 - 1, (line.Length >= 103 ? ((103 - 99) + 1) : line.Length - (99 - 1))).Trim() : "";
                PHI = (105 > 0 && 109 > 0 && line.Length >= 105) ? line.Substring(105 - 1, (line.Length >= 109 ? ((109 - 105) + 1) : line.Length - (105 - 1))).Trim() : "";
                PSI = (111 > 0 && 115 > 0 && line.Length >= 111) ? line.Substring(111 - 1, (line.Length >= 115 ? ((115 - 111) + 1) : line.Length - (111 - 1))).Trim() : "";
                X_CA = (116 > 0 && 122 > 0 && line.Length >= 116) ? line.Substring(116 - 1, (line.Length >= 122 ? ((122 - 116) + 1) : line.Length - (116 - 1))).Trim() : "";
                Y_CA = (123 > 0 && 129 > 0 && line.Length >= 123) ? line.Substring(123 - 1, (line.Length >= 129 ? ((129 - 123) + 1) : line.Length - (123 - 1))).Trim() : "";
                Z_CA = (130 > 0 && 136 > 0 && line.Length >= 130) ? line.Substring(130 - 1, (line.Length >= 136 ? ((136 - 130) + 1) : line.Length - (130 - 1))).Trim() : "";

                if (char.IsWhiteSpace(SecondaryStructure)) SecondaryStructure = dssp_default_secondary_structure;
            }

        }

        public static List<DsspRecord> Load(string file)
        {
            const string SsDataMarker1 = "  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA ";
            const string SsDataMarker2 = "  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA            CHAIN";
            var lines = File.ReadAllLines(file);
            
            var dssp_record_list = new List<DsspRecord>();

            var marker_found = false;

            foreach (var line in lines)
            {
                if (!marker_found && (line == SsDataMarker1 || line == SsDataMarker2))
                {
                    marker_found = true;
                    continue;
                }

                if (marker_found)
                {
                    var dr = new DsspRecord(line);
                    dssp_record_list.Add(dr);
                }
            }

            return dssp_record_list;
        }
    }
}