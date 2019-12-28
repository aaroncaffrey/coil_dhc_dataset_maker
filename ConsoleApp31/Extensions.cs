using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Coil_DHC_Dataset_Maker
{
    public static class Extensions
    {

        public static T[][] ToJagged<T>(this T[,] value)
        {
            if (Object.ReferenceEquals(null, value))
                return null;

            // Jagged array creation
            T[][] result = new T[value.GetLength(0)][];

            for (int i = 0; i < value.GetLength(0); ++i)
                result[i] = new T[value.GetLength(1)];

            // Jagged array filling
            for (int i = 0; i < value.GetLength(0); ++i)
            for (int j = 0; j < value.GetLength(1); ++j)
                result[i][j] = value[i, j];

            return result;
        }
        public static int? NullableTryParseInt32(this string str)
        {
            int intValue;
            return int.TryParse(str, out intValue) ? (int?)intValue : null;
        }

        public static double try_parse_double(this string value, double default_value = 0)
        {
            if (!string.IsNullOrWhiteSpace(value) && double.TryParse(value, NumberStyles.Float, CultureInfo.InvariantCulture, out var result))
            {
                return result;
            }
            else
            {
                return default_value;
            }
        }

        public static List<int> IndexOfAll(this string text, char[] substring)
        {
            return substring.SelectMany(a => IndexOfAll(text, a)).ToList();
        }

        public static List<int> IndexOfAll(this string text, List<char> substring)
        {
            return substring.SelectMany(a => IndexOfAll(text, a)).ToList();
        }

        public static List<int> IndexOfAll(this string text, char substring)
        {
            var indexes = new List<int>();

            if (string.IsNullOrEmpty(text)) return new List<int>();
            
            int index = -1;

            do
            {
                index = text.IndexOf(substring, index + 1);//, StringComparison.InvariantCulture);
                if (index != -1) indexes.Add(index);

            } while (index != -1);
        

            return indexes;
        }

        public static List<int> IndexOfAll(this string text, string substring)
        {
            var indexes = new List<int>();

            if (string.IsNullOrEmpty(text)) return new List<int>();


            if (!string.IsNullOrEmpty(text) && !string.IsNullOrEmpty(substring))
            {
                int index = -1;

                do
                {
                    index = text.IndexOf(substring, index + 1, StringComparison.InvariantCulture);
                    if (index != -1) indexes.Add(index);

                } while (index != -1);
            }

            return indexes;
        }
    }
}
