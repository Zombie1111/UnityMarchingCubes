#if UNITY_EDITOR
using UnityEngine;
using UnityEditor;
using System.IO;
using System.Text.RegularExpressions;

namespace zombGen_editor
{
    //AI generated, used to help generating GetTriangleConnectionTables array
    public class FindAndIncreasement : EditorWindow
    {
        private string patternText = "";
        private int increaseBy = 1;
        private int startValue = 0;

        [MenuItem("Tools/Helpers/Find And Increasement")]
        public static void ShowWindow()
        {
            GetWindow<FindAndIncreasement>("Find And Increasement");
        }

        private void OnGUI()
        {
            GUILayout.Label("Find And Increasement Tool", EditorStyles.boldLabel);
            EditorGUILayout.Space();

            patternText = EditorGUILayout.TextField("Text to Match", patternText);
            increaseBy = EditorGUILayout.IntField("Increase By", increaseBy);
            startValue = EditorGUILayout.IntField("Start Value", startValue);

            EditorGUILayout.Space();

            if (GUILayout.Button("Execute"))
            {
                ExecuteOperation();
            }
        }

        private void ExecuteOperation()
        {
            if (string.IsNullOrEmpty(patternText))
            {
                EditorUtility.DisplayDialog("Error", "Please enter a match text containing a number.", "OK");
                return;
            }

            string filePath = EditorUtility.OpenFilePanel("Select Text File", "", "txt");
            if (string.IsNullOrEmpty(filePath))
                return;

            string fileContent = File.ReadAllText(filePath);

            // Escape user input for regex, but allow digits placeholder
            string escapedPattern = Regex.Escape(patternText);
            // Find the first number sequence in the text to replace dynamically
            escapedPattern = Regex.Replace(escapedPattern, @"\\d+", @"(\d+)");
            Regex regex = new Regex(escapedPattern);

            int count = 0;
            string result = regex.Replace(fileContent, match =>
            {
                int newValue = startValue + (count * increaseBy);
                count++;
                string original = match.Value;
                // Replace only the numeric part
                return Regex.Replace(original, @"\d+", newValue.ToString());
            });

            File.WriteAllText(filePath, result);
            EditorUtility.DisplayDialog("Success", $"Updated {count} matches and saved back to file.", "OK");
        }
    }
}
#endif
