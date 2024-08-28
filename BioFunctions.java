import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Scanner;

public class BioFunctions {
    // Methods to reuse in multiple Code Challenges
    // (input/output, helper methods etc - not the actual work)

    public static char[] nucleotides = {'A', 'C', 'G', 'T'};

    // Method to get lines from a file into an ArrayList of Strings
    public static ArrayList<String> readFile(String filename) {
        try {
            Scanner reader = new Scanner(new File(filename));
            ArrayList<String> list = new ArrayList<>();

            while (reader.hasNextLine()) {
                list.add(reader.nextLine());
            }
            reader.close();

            return list;

        } catch (FileNotFoundException e) {
            System.err.println(e.getMessage());
            return null;
        }
    }

    // Helper method that returns the base pair of a DNA nucleotide (A-T, C-G)
    // Returns '-' if invalid letter
    public static char basePair(char nucleotide) {
        switch (nucleotide) {
            case ('A'):
                return 'T';
            case ('T'):
                return 'A';
            case ('C'):
                return 'G';
            case ('G'):
                return 'C';
            default:
                return '-';
        }
    }

     public static int baseNum(char letter) {
        switch (letter) {
            case ('A'):
                return 0;

            case ('C'):
                return 1;

            case ('G'):
                return 2;

            case ('T'):
                return 3;
        }
        return -1;
    }

    // Hamming distance: number of 'mismatches' between strings
    // for example abcd and abed have a hd of 1
    public static int hammingDistance(String one, String two) throws IllegalArgumentException {
        if (one.length() != two.length()) {
            throw new IllegalArgumentException("Strings must be the same length");
        }
        int count = 0;
        for (int i = 0; i < one.length(); i++) {
            if (one.charAt(i) != two.charAt(i)) {
                count++;
            }
        }
        return count;
    }

    // literally just a printing method 
    // works with Lists and Sets
    public static void print(Collection<String> text) {
        for (Object word : text) {
            System.out.print(word + " ");
        }
        System.out.println();
    }

   
}
