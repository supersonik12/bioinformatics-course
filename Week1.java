import java.util.HashMap;
import java.util.HashSet;
import java.util.ArrayList;
import java.util.Set;

public class Week1 {
    // Coding Challenges from Week 1 of Bioinformatics I

    public static void main(String[] args) {
        ArrayList<String> file = BioFunctions.readFile("Week 1 Datasets/dataset_4_5.txt");
        
        System.out.println(findClumps(file.get(0), 9, 500, 3).size());

        System.out.println(patternMatch("ATA", "GACGATATACGACGATA"));
    }

    

    // Counts instances of a String pattern in String text (can overlap) 
    // Finding repeating k-mers in regions of DNA to locate ORI
    public static int patternCount(String text, String pattern) {
        int count = 0;

        for (int i = 0; i <= text.length()-pattern.length(); i++) {
            if (text.substring(i, i + pattern.length()).equals(pattern)) {
                count++;
            }
        }
        return count;

    }

    // Finds most frequent strings of length k in text
    // Finding the most frequent k-mers in a region of DNA given the length
    public static ArrayList<String> frequentWords(String text, int k) {
        ArrayList<String> frequent = new ArrayList<>();
        HashMap<String, Integer> map = frequencyTable(text, k);
        
        int max = (int) map.values().toArray()[0];
        for (Integer count : map.values()) {
            if (count > max) {
                max = count;
            }
        }

        for (String pattern : map.keySet()) {
            if (map.get(pattern) == max) {
                frequent.add(pattern);
            }
        }

        return frequent;
    }

    // Gets the frequencies of every string of length k in text
    // (used to be part of frequentWords until I needed it for findClumps)
    public static HashMap<String, Integer> frequencyTable(String text, int k) {
        HashMap<String, Integer> map = new HashMap<>();
        // add values to map
        for (int i = 0; i <= text.length()-k; i++) {
            String pattern = text.substring(i, i+k);
            if (map.containsKey(pattern)) {
                map.put(pattern, map.get(pattern) + 1);
            } else {
                map.put(pattern, 1);
            }
        }
        return map;
    }

    // Finds the complement String of DNA using base pairs
    // Returns string in the 5' to 3' direction
    public static String reverseComplement(String input) {
        String output = "";

        for (int i = input.length()-1; i >= 0; i--) {
            output += BioFunctions.basePair(input.charAt(i));
        }

        return output;
    }
    
    // Finds the starting postitions of every spot pattern exists in genome
    public static ArrayList<Integer> patternMatch(String pattern, String genome) {
        ArrayList<Integer> output = new ArrayList<>();
        for (int i = 0; i <= genome.length()-pattern.length(); i++) {
            if (genome.substring(i, i+pattern.length()).equals(pattern)) {
                output.add(i);
            }
        }
        return output;
    }

    // Finds strings of length k in genome that repeat t times within L characters
    // Finds spots where repeating k-mers exist in the genome (clumps) that indicate ORI
    public static Set<String> findClumps(String genome, int k, int L, int t) {
        // Fairly inefficient but works for now
        // better way - iterate once, count every instance of strings of length k
        // then go over that, check if instances are within L of each other and > t

        // Using a set to avoid duplicate strings
        Set<String> set = new HashSet<>();
        
        for (int i = 0; i <= genome.length()-L; i++) {
            String window = genome.substring(i, i+L);
            HashMap<String, Integer> map = frequencyTable(window, k);
            for (String pattern : map.keySet()) {
                if (map.get(pattern) >= t) {
                    set.add(pattern);
                }
            }
        }

        return set;

    }
}