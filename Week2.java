import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class Week2 {
    public static void main(String[] args) {
        ArrayList<String> input = BioFunctions.readFile("Week 2 Datasets/Salmonella_enterica.txt");
        String genome = input.get(1);
        System.out.println("Minimum skew at positions " + minSkew(genome));
        System.out.println("Predicted ORI +- 500 bases: " + genome.substring(3764856-500, 3764856+500));
        System.out.println("Most frequent 9-mers (within HD of 1):\n" + approxFreqWords(genome.substring(3764856-500, 3764856+500), 9, 1));
        System.out.println(approxMatch("TTATCCACA", genome.substring(3764856-500, 3764856+500), 1));

        // Most frequent 9-mers (within HD of 1): [TTATCCACA, TGTGGATAA]
        // Looks like S. enterica has the same DnaA box as E. coli - TTATCCACA
    }



    // Counts difference between 'C' and 'G' in genome
    public static int[] skew(String genome, int i) {
        int[] skews = new int[i+1];
        int c = 0;
        int g = 0;
        for (int index = 0; index < i; index++) {
            // Start at 0, add one every time
            //System.out.print(g - c + " ");
            skews[index] = g - c;
            if (genome.charAt(index) == 'C') {
                c++;
                // no, java
            } else if (genome.charAt(index) == 'G') {
                g++;
            }
        }
        //System.out.println(g - c);
        skews[i] = g - c;

        return skews;
    }

    // Finding all indices of minimum skew in genome
    public static ArrayList<Integer> minSkew(String genome) {
        int[] skews = skew(genome, genome.length());
        ArrayList<Integer> mins = new ArrayList<>();
        int min = skews[0];
        for (int i = 0; i < skews.length; i++) {
            if (skews[i] == min) {
                min = skews[i];
                mins.add(i);
            } else if (skews[i] < min) {
                mins.clear();
                min = skews[i];
                mins.add(i);
            }
        }
        return mins;
    }

    
    // Find approximate matches to strings within a Hamming distance of int hd
    // includes reverse complements (as negative numbers for now)
    public static ArrayList<Integer> approxMatch(String pattern, String genome, int hd) {
         ArrayList<Integer> output = new ArrayList<>();
        for (int i = 0; i <= genome.length()-pattern.length(); i++) {
            String current = genome.substring(i, i+pattern.length());
            if (BioFunctions.hammingDistance(current, pattern) <= hd) {
                output.add(i);
                //System.out.print(i + " ");
            } 
            else if (BioFunctions.hammingDistance(current, Week1.reverseComplement(pattern)) <= hd) {
                output.add(-i);
                //System.out.print(i + " ");
            }
        }
        return output;
    }

    // i dont know how this works but it does
    // figure out later
    public static HashSet<String> neighbors(String pattern, int distance) {
        // this is going to be recursive i can just feel it
        HashSet<String> neighborhood = new HashSet<>();
        if (distance == 0) {
            neighborhood.add(pattern);
            return neighborhood;
        }
        if (pattern.length() == 1) {
            neighborhood.add("A");
            neighborhood.add("C");
            neighborhood.add("G");
            neighborhood.add("T");
            return neighborhood;
        }

        HashSet<String> suffixNeighbors = neighbors(pattern.substring(1), distance);

        for (String suffixed : suffixNeighbors) {
            if (BioFunctions.hammingDistance(suffixed, pattern.substring(1)) < distance) {
                for (char nt : BioFunctions.nucleotides) {
                    neighborhood.add(nt + suffixed);
                }
            } else {
                neighborhood.add(pattern.charAt(0) + suffixed);
            }
        } 
        return neighborhood;
    }

    // Find the most frequent k-mers within a Hamming distance of d
    // FrequentWords from Week1 with mismatches and reverse complements
    public static ArrayList<String> approxFreqWords(String text, int k, int d) {
        ArrayList<String> frequent = new ArrayList<>();
        HashMap<String, Integer> map =  new HashMap<>();

        // add values to frequency map
        for (int i = 0; i <= text.length()-k; i++) {
            String pattern = text.substring(i, i+k);
            HashSet<String> neighborhood = neighbors(pattern, d);
            for (String neighbor : neighborhood) {
                if (map.containsKey(neighbor)) {
                    map.put(neighbor, map.get(neighbor) + 1);
                } else {
                    map.put(neighbor, 1);
                }
                // reverse complements here
                String reversed = Week1.reverseComplement(neighbor);
                if (map.containsKey(reversed)) {
                    map.put(reversed, map.get(reversed) + 1);
                } else {
                    map.put(reversed, 1);
                }
            }
        }

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

}
