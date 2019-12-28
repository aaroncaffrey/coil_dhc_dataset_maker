using System.Collections.Generic;
using System.IO;

namespace Coil_DHC_Dataset_Maker
{
    public static class ColumnValueExtensions
    {
        public static string ColumnValue(this string columnFormatLine, int startPosition, int endPosition, bool trim = true)
        {
            var r = columnFormatLine.Substring(startPosition - 1, (endPosition - startPosition) + 1);
            if (trim) r = r.Trim();
            return r;
        }
    }

    public class StrideFile
    {
        public abstract class Stride_Record
        {
            //1-3	Record code
            //4-5	Not used
            //6-73	Data
            //74-75	Not used
            //75-79	Four letter PDB code (if available)

            public string RecordCode;
            public string NotUsed1;
            public string Data;
            public string NotUsed2;
            public string PdbCode;

            protected Stride_Record(string columnFormatLine)
            {
                if (string.IsNullOrWhiteSpace(columnFormatLine)) return;

                RecordCode = columnFormatLine.ColumnValue(1, 3);
                NotUsed1 = columnFormatLine.ColumnValue(4, 5);
                Data = columnFormatLine.ColumnValue(6, 73);
                NotUsed2 = columnFormatLine.ColumnValue(74, 75);
                PdbCode = columnFormatLine.ColumnValue(75, 79);
            }
        }

        public class Stride_Remark : Stride_Record
        {
            //REM Remarks and blank   lines

            public Stride_Remark(string columnFormatLine) : base(columnFormatLine)
            {

            }
        }

        public class Stride_Header : Stride_Record
        {
            //HDR Header.  Protein name, date of file creation and PDB code

            public Stride_Header(string columnFormatLine) : base(columnFormatLine)
            {

            }
        }

        public class Stride_Compound : Stride_Record
        {
            //CMP Compound.Full name of the molecule and identifying information

            public Stride_Compound(string columnFormatLine) : base(columnFormatLine)
            {

            }
        }

        public class Stride_Receptor : Stride_Record
        {
            //SRC Species, organ, tissue, and mutant from which the molecule has been obtained

            public Stride_Receptor(string columnFormatLine) : base(columnFormatLine)
            {

            }
        }

        public class Stride_StructureAuthors : Stride_Record
        {
            //AUT Names of the structure authors

            public Stride_StructureAuthors(string columnFormatLine) : base(columnFormatLine)
            {

            }

        }

        public class Stride_Chain : Stride_Record
        {
            //CHN File name and PDB chain identifier*).
            // Format: File name beginning from position 6 followed by one space and one-letter chain    identifier

            // 0        1         2         3         4         5         6         7
            // 1234567890123456789012345678901234567890123456789012345678901234567890123456789
            // CHN  /webclu/data/stride/pdb/pdb1a12.ent C                                 1A12

            public string FileName;
            public string OneLetterChainIdentifier;

            public Stride_Chain(string columnFormatLine) : base(columnFormatLine)
            {
                FileName = columnFormatLine.ColumnValue(6, 75);
                OneLetterChainIdentifier = FileName.Substring(FileName.LastIndexOf(' ') + 1);
                FileName = FileName.Substring(0, FileName.LastIndexOf(' '));
            }
        }

        public class Stride_AminoAcidSequence : Stride_Record
        {
            //SEQ Amino acid sequence
            //
            //Format:
            //  6-9  First residue PDB number
            // 11-60 Sequence
            // 62-65 Last residue PDB number

            public string FirstResiduePdbNumber;
            public string Sequence;
            public string LastResiduePdbNumber;

            public Stride_AminoAcidSequence(string columnFormatLine) : base(columnFormatLine)
            {
                FirstResiduePdbNumber = columnFormatLine.ColumnValue(6, 9);
                Sequence = columnFormatLine.ColumnValue(11, 60);
                LastResiduePdbNumber = columnFormatLine.ColumnValue(62, 65);
            }
        }

        public class Stride_SecondaryStructure : Stride_Record
        {

            //STR Secondary   structure summary
            //
            //Format:
            // 11-60 Secondary structure assignment**

            public string SecondaryStructureAssignment;

            public Stride_SecondaryStructure(string columnFormatLine) : base(columnFormatLine)
            {
                SecondaryStructureAssignment = columnFormatLine.ColumnValue(11, 60);
            }
        }

        public class Stride_LocationOfSecondaryStructureElements : Stride_Record
        {
            //LOC Location of secondary structure elements
            //
            //Format:
            //  6-17 Element name
            // 19-21 First residue name
            // 22-27 First residue PDB number (spec says 32-26 !)
            // 28-28 First residue chain identifier
            // 36-38 Last residue name
            // 42-45 Last residue PDB number
            // 47-47 Last residue chain identifier

            // 0        1         2         3         4         5         6         7
            // 1234567890123456789012345678901234567890123456789012345678901234567890123456789
            // LOC  GammaInv     LEU   117 C      GLU    119 C                            1A12
            public string ElementName;
            public string FirstResidueName;
            public string FirstResiduePdbNumber;
            public string FirstResidueChainIdentifier;
            public string LastResidueName;
            public string LastResiduePdbNumber;
            public string LastResidueChainIdentifier;

            public Stride_LocationOfSecondaryStructureElements(string columnFormatLine) : base(columnFormatLine)
            {
                ElementName = columnFormatLine.ColumnValue(6, 17);
                FirstResidueName = columnFormatLine.ColumnValue(19, 21);
                FirstResiduePdbNumber = columnFormatLine.ColumnValue(22, 27);
                FirstResidueChainIdentifier = columnFormatLine.ColumnValue(28, 28);
                LastResidueName = columnFormatLine.ColumnValue(36, 38);
                LastResiduePdbNumber = columnFormatLine.ColumnValue(42, 45);
                LastResidueChainIdentifier = columnFormatLine.ColumnValue(47, 47);
            }
        }

        public class Stride_DetailedSecondaryStructureAssignments : Stride_Record
        {
            //ASG Detailed secondary structure assignment
            //
            //Format:
            //  6-8  Residue name
            // 10-10 Protein chain identifier
            // 12-15 PDB residue number
            // 17-20 Ordinal residue number
            // 25-25 One letter secondary structure code**)
            // 27-39 Full secondary structure name
            // 43-49 Phi angle
            // 53-59 Psi angle
            // 65-69 Residue solvent accessible area

            public const char stride_default_secondary_structure = 'C';

            public string ResidueName;
            public string ProteinChainIdentifier;
            public string PdbResidueNumber;
            public string OrdinalResidueNumber;
            public string OneLetterSecondaryStructureCode;
            public string FullSecondaryStructureName;
            public string PhiAngle;
            public string PsiAngle;
            public string ResidueSolventAccessibleArea;

            public Stride_DetailedSecondaryStructureAssignments(string columnFormatLine) : base(columnFormatLine)
            {
                ResidueName = columnFormatLine.ColumnValue(6, 8);
                ProteinChainIdentifier = columnFormatLine.ColumnValue(10, 10);
                PdbResidueNumber = columnFormatLine.ColumnValue(12, 15);
                OrdinalResidueNumber = columnFormatLine.ColumnValue(17, 20);
                OneLetterSecondaryStructureCode = columnFormatLine.ColumnValue(25, 25);
                FullSecondaryStructureName = columnFormatLine.ColumnValue(27, 39);
                PhiAngle = columnFormatLine.ColumnValue(43, 49);
                PsiAngle = columnFormatLine.ColumnValue(53, 59);
                ResidueSolventAccessibleArea = columnFormatLine.ColumnValue(65, 69);

                if (string.IsNullOrWhiteSpace(OneLetterSecondaryStructureCode)) OneLetterSecondaryStructureCode = stride_default_secondary_structure.ToString();
            }

        }

        public class Stride_LigandDonorReside : Stride_Record
        {
            // DNR LigandDonor residue
            //
            //Format:
            // 6-8  LigandDonor residue name
            //10-10 Protein chain identifier
            //12-15 PDB residue number
            //17-20 Ordinal residue number
            //26-28 Acceptor residue name
            //30-30 Protein chain identifier
            //32-35 PDB residue number
            //37-40 Ordinal residue number
            //42-45 N__0 distance
            //47-52 N__O_C angle
            //54-59 O__N_C angle
            //61-66 Angle between the planes of ligand_donor complex and O__N_C
            //68-73 Angle between the planes of acceptor complex and N__O_C

            public string LigandDonorResidueName;
            public string ProteinChainIdentifier1;
            public string PdbResidueNumber1;
            public string OrdinalResidueNumber1;
            public string AcceptorResidueName;
            public string ProteinChainIdentifier2;
            public string PdbResidueNumber2;
            public string OrdinalResidueNumber2;
            public string N__0Distance;
            public string N__o_cAngle;
            public string O__n_cAngle;
            public string AngleBetweenThePlanesOfLigandDonorComplexAndO__n_c;
            public string AngleBetweenThePlanesOfAcceptorComplexAndN__o_c;

            public Stride_LigandDonorReside(string columnFormatLine) : base(columnFormatLine)
            {
                LigandDonorResidueName = columnFormatLine.ColumnValue(6, 8);
                ProteinChainIdentifier1 = columnFormatLine.ColumnValue(10, 10);
                PdbResidueNumber1 = columnFormatLine.ColumnValue(12, 15);
                OrdinalResidueNumber1 = columnFormatLine.ColumnValue(17, 20);
                AcceptorResidueName = columnFormatLine.ColumnValue(26, 28);
                ProteinChainIdentifier2 = columnFormatLine.ColumnValue(30, 30);
                PdbResidueNumber2 = columnFormatLine.ColumnValue(32, 35);
                OrdinalResidueNumber2 = columnFormatLine.ColumnValue(37, 40);
                N__0Distance = columnFormatLine.ColumnValue(42, 45);
                N__o_cAngle = columnFormatLine.ColumnValue(47, 52);
                O__n_cAngle = columnFormatLine.ColumnValue(54, 59);
                AngleBetweenThePlanesOfLigandDonorComplexAndO__n_c = columnFormatLine.ColumnValue(61, 66);
                AngleBetweenThePlanesOfAcceptorComplexAndN__o_c = columnFormatLine.ColumnValue(68, 73);
            }
        }

        public class Stride_AcceptorResidue : Stride_Record
        {
            // ACC Acceptor residue
            //
            //Format:
            // 6-8  Acceptor residue name
            //10-10 Protein chain identifier
            //12-15 PDB residue number
            //17-20 Ordinal residue number
            //26-28 LigandDonor residue name
            //30-30 Protein chain identifier
            //32-35 PDB residue number
            //37-40 Ordinal residue number
            //42-45 N..0 distance
            //47-52 N..O=C angle
            //54-59 O..N-C angle
            //61-66 Angle between the planes of ligand_donor complex and O..N-C
            //68-73 angle between the planes of acceptor complex and N..O=C

            public string AcceptorResidueName;
            public string ProteinChainIdentifier1;
            public string PdbResidueNumber1;
            public string OrdinalResidueNumber1;
            public string LigandDonorResidueName;
            public string ProteinChainIdentifier2;
            public string PdbResidueNumber2;
            public string OrdinalResidueNumber2;
            public string N__0Distance;
            public string N__o_cAngle;
            public string O__n_cAngle;
            public string AngleBetweenThePlanesOfLigandDonorComplexAndO__n_c;
            public string AngleBetweenThePlanesOfAcceptorComplexAndN__o_c;

            public Stride_AcceptorResidue(string columnFormatLine) : base(columnFormatLine)
            {
                AcceptorResidueName = columnFormatLine.ColumnValue(6, 8);
                ProteinChainIdentifier1 = columnFormatLine.ColumnValue(10, 10);
                PdbResidueNumber1 = columnFormatLine.ColumnValue(12, 15);
                OrdinalResidueNumber1 = columnFormatLine.ColumnValue(17, 20);
                LigandDonorResidueName = columnFormatLine.ColumnValue(26, 28);
                ProteinChainIdentifier2 = columnFormatLine.ColumnValue(30, 30);
                PdbResidueNumber2 = columnFormatLine.ColumnValue(32, 35);
                OrdinalResidueNumber2 = columnFormatLine.ColumnValue(37, 40);
                N__0Distance = columnFormatLine.ColumnValue(42, 45);
                N__o_cAngle = columnFormatLine.ColumnValue(47, 52);
                O__n_cAngle = columnFormatLine.ColumnValue(54, 59);
                AngleBetweenThePlanesOfLigandDonorComplexAndO__n_c = columnFormatLine.ColumnValue(61, 66);
                AngleBetweenThePlanesOfAcceptorComplexAndN__o_c = columnFormatLine.ColumnValue(68, 73);
            }
        }


        public static List<Stride_Record> Load(string strideFilename)
        {
            var result = new List<Stride_Record>();

            var lines = File.ReadAllLines(strideFilename);

            foreach (var columnFormatLine in lines)
            {
                if (string.IsNullOrWhiteSpace(columnFormatLine) || columnFormatLine.Length <= 3) continue;

                var recordCode = columnFormatLine.Substring(0, 3);
                Stride_Record record = null;

                switch (recordCode)
                {
                    case "REM":
                        record = new Stride_Remark(columnFormatLine); break;
                    case "HDR":
                        record = new Stride_Header(columnFormatLine); break;
                    case "CMP":
                        record = new Stride_Compound(columnFormatLine); break;
                    case "SRC":
                        record = new Stride_Receptor(columnFormatLine); break;
                    case "AUT":
                        record = new Stride_StructureAuthors(columnFormatLine); break;
                    case "CHN":
                        record = new Stride_Chain(columnFormatLine); break;
                    case "SEQ":
                        record = new Stride_AminoAcidSequence(columnFormatLine); break;
                    case "STR":
                        record = new Stride_SecondaryStructure(columnFormatLine); break;
                    case "LOC":
                        record = new Stride_LocationOfSecondaryStructureElements(columnFormatLine); break;
                    case "ASG":
                        record = new Stride_DetailedSecondaryStructureAssignments(columnFormatLine); break;
                    case "DNR":
                        record = new Stride_LigandDonorReside(columnFormatLine); break;
                    case "ACC":
                        record = new Stride_AcceptorResidue(columnFormatLine); break;
                    default:
                        //record = new Stride_Record(columnFormatLine);
                        break;
                }

                if (record != null) result.Add(record);
            }

            return result;
        }

    }
}
