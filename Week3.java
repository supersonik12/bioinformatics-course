import java.util.HashSet;
import java.util.ArrayList;
import java.util.Arrays;

public class Week3 {
    public static void main(String[] args) {
        ArrayList<String> input = BioFunctions.readFile("Week 3 Datasets/subtle_motif_dataset.txt");
  
        // Solving the Subtle Motif Problem
        // -> Greedy Motif Search implmenting Laplace's Rule  
        String[] DNA = input.toArray(new String[0]);
        String[] motifs = greedyMotifSearch(DNA, 15);
        System.out.println(Arrays.toString(motifs));
        System.out.println("Implanted Motif was: AAAAAAAAGGGGGGG");
        System.out.println("Consensus String:    " + bestMotif(motifs));
        System.out.println("Score: " + score(motifs));

    }

    // very very inefficient
    // I felt bad writing this
    public static HashSet<String> motifEnumeration(String[] DNA, int k, int d) {
        HashSet<String> patterns = new HashSet<>();

        for (int i = 0; i <= DNA[0].length() - k; i++) {
            HashSet<String> neighbors = Week2.neighbors(DNA[0].substring(i, i + k), d);
            for (String neighbor : neighbors) {
                boolean contains = true;
                for (String gene : DNA) {
                    if (!checkForMatch(neighbor, gene, d)) {
                        contains = false;
                        break;
                    }
                }
                if (contains) {
                    patterns.add(neighbor);
                }
            }
        }
        return patterns;
    }

    public static boolean checkForMatch(String pattern, String genome, int hd) {
        for (int i = 0; i <= genome.length() - pattern.length(); i++) {
            String current = genome.substring(i, i + pattern.length());
            if (BioFunctions.hammingDistance(current, pattern) <= hd) {
                return true;
            }

        }
        return false;
    }

    // --- Finding an 'ideal' motif from a matrix of motifs ---
    public static String bestMotif(String[] motifs) {
        int[][] count = count(motifs);
        double[][] profile = profile(count, motifs.length);
        return consensus(profile);
    }

    // Count occurences of a base and stores in an array
    // Array will always be 4 x motifArray length
    public static int[][] count(String[] motifArray) {
        int[][] countArray = new int[4][motifArray[0].length()];

        for (int col = 0; col < motifArray[0].length(); col++) {
            for (int row = 0; row < motifArray.length; row++) {
                int letter = BioFunctions.baseNum(motifArray[row].charAt(col));
                countArray[letter][col] += 1;
                
            }

        }
        return countArray;
    }

    // Creates an array representing how likely a specific base is
    // 'Probability' of a certain base occuring at that index in motif
    // (version implementing pseudocounts (Laplace Rule) in next section)
    public static double[][] profile(int[][] countArray, int numOfMotifs) {
        double[][] profile = new double[4][countArray[0].length];

        for (int row = 0; row < countArray.length; row++) {
            for (int col = 0; col < countArray[row].length; col++) {
                profile[row][col] = (double) countArray[row][col] / numOfMotifs;
            }
        }

        return profile;
    }

    public static String consensus(double[][] profile) {
        String output = "";
        for (int col = 0; col < profile[0].length; col++) {
            double max = 0;
            int maxIndex = 0;
            for (int row = 0; row < profile.length; row++) {
                if (profile[row][col] > max) {
                    max = profile[row][col];
                    maxIndex = row;
                }
            }
            output += BioFunctions.nucleotides[maxIndex];
        }
        return output;
    }

    // Scores an array based on Hamming Distance
    // Finds how many occurences of NOT the most common base
    // Lower scores = arrays are closer to each other
    public static int score(String[] motifArray) {
        int[][] countArray = count(motifArray);
        int score = 0;

        for (int i = 0; i < countArray[0].length; i++) {
            int a = Math.max(countArray[0][i], countArray[1][i]);
            int b = Math.max(countArray[2][i], countArray[3][i]);
            score += motifArray.length - Math.max(a, b);
        }
        return score;
    }



    // --- Motif finding using the 'median string' method

    /*MedianString(Dna, k)
    distance ← ∞
    for each k-mer Pattern from AA…AA to TT…TT
        if distance > d(Pattern, Dna)
             distance ← d(Pattern, Dna)
             Median ← Pattern
    return Median */
    /* for each string Text in Dna
        HammingDistance ← ∞
        for each k-mer Pattern’ in Text
            if HammingDistance > HammingDistance(Pattern, Pattern’)
                HammingDistance ← HammingDistance(Pattern, Pattern’)
        distance ← distance + HammingDistance */

    public static int hammingSum(String pattern, String[] motifs) {
        int k = pattern.length();
        int sum = 0;
        for (String motif : motifs) {
            int distance = pattern.length();
            for (int i = 0; i <= motif.length()- k; i++) {
                int d = BioFunctions.hammingDistance(pattern, motif.substring(i, i+k));
                if (d < distance) {
                    distance = d;
                }
            }
            sum += distance;
        }
        return sum;
    }

    public static String medianString(String[] DNA, int k) {
        int distance = k * DNA.length;
        String median = "";

        // Generate all possible strings of length k 
        // { A C G T }, uses neighbors method
        String s = "";
        for (int i = 0; i < k; i++) {
            s += "A";
        }
        HashSet<String> kmers = Week2.neighbors(s, k);

        for (String pattern : kmers) {
            int d = hammingSum(pattern, DNA);
            if (d < distance) {
                distance = d;
                median = pattern;
            }
        }   
        return median;     
    }

    // --- Greedy motif search using probable strings ---4

    // TODO: more comments
    public static String[] greedyMotifSearch(String[] DNA, int k) {
        int t = DNA.length;
       
        String[] bestMotifs = new String[t];
        for (int count = 0; count < t; count++) {
            bestMotifs[count] = DNA[count].substring(0, k);
        }
        // position - of a kmer in a string of DNA
        // index - of string in DNA array 
        for (int position = 0; position <= DNA[0].length()-k; position++) {
            ArrayList<String> currentMotifList = new ArrayList<>();
            currentMotifList.add(DNA[0].substring(position, position + k));
            for (int index = 1; index < t; index++) {
                // 'Worse' version used profile(), which messes up when 0 is involved
                double[][] profile = laplaceProfile(count((String[]) currentMotifList.toArray(new String[0])),
                                                    currentMotifList.size());
                currentMotifList.add(mostProbableString(DNA[index], k, profile));
            }
            if (score((String[]) currentMotifList.toArray(new String[0])) < score(bestMotifs)) {
                bestMotifs = (String[]) currentMotifList.toArray(new String[0]);
            }
        }
        return bestMotifs;
    }

    // Calculate the probability of a motif occuring according to a given profile
    // FOR BETTER RESULTS INCLUDE PSEUDOCOUNTS (see Laplace's Rule of Succession)
    public static double probabilty(String motif, double[][] profile) {
        if (profile[0].length != motif.length()) {
            throw new IllegalArgumentException("Motif length doesn't match profile length");
        }
        double probability = 1;
        for (int i = 0; i < motif.length(); i++) {
            int number = BioFunctions.baseNum(motif.charAt(i));
            probability *= profile[number][i];

        }
        return probability;
    }

    public static String mostProbableString(String text, int k, double[][] profile) {
        double highest = 0.0;
        String output = text.substring(0, k);
        for (int i = 0; i <= text.length() - k; i++) {
            String current = text.substring(i, i+k);
            double probability = probabilty(current, profile);
            if (probability > highest) {
                output = current;
                highest = probability;
            }
        }
        return output;
    }

    public static String[] mostProbableArray(String[] DNA, int k, double[][] profile) {
        String[] found = new String[DNA.length];
        for (int i = 0; i < DNA.length; i++) {
            found[i] = mostProbableString(DNA[i], k, profile);
        }
        return found;
    }

    // Works just like profile() but implements pseudocounts to account for 
    // 'rare but not impossible' bases
    // Prevents the 0-cancels-everything-out effect
    public static double[][] laplaceProfile(int[][] countArray, int numOfMotifs) {
        double[][] profile = new double[4][countArray[0].length];

        for (int row = 0; row < countArray.length; row++) {
            for (int col = 0; col < countArray[row].length; col++) {
                countArray[row][col] += 1;
                profile[row][col] = (double) countArray[row][col] / (numOfMotifs + 4);
            }
        }

        return profile;
    }
}
