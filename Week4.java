import java.util.Arrays;

public class Week4 {

    private static boolean tracing = false;
    // Randomized Approaches to Finding Motifs

    public static void main(String[] args) {
        // For coding problems:
        // ArrayList<String> input = BioFunctions.readFile("Week 4 Datasets/dataset_163_4.txt");
        // String[] DNA = input.get(1).split(" ");

        // Subtle Motif Problem:
        String[] DNA = BioFunctions.readFile("Week 4 Datasets/DosR_dataset.txt").toArray(new String[0]);
        
        // --- Randomized Motif Search ---
        long start = System.currentTimeMillis();
        String[] motifs = randomizedMotifSearch(DNA, 15, 10000);
        long end = System.currentTimeMillis();
        
        System.out.println("Best match: " + Arrays.toString(motifs));
        System.out.println("Consensus String: " + Week3.bestMotif(motifs));
        System.out.println("Score: " + Week3.score(motifs));
        System.out.println((end - start)/1000 + " s");

        // ran w 100,000 trials, came up with this -
        // Best: {"CAGAAACGAGAGGAG", "TAAAAAATAGCAGGG", "AAATAAACAGCGGGG", "ACAGAAAAAAAGGGG", "CATAAAGTAGAGGGG", "CTAAAAATGGGGCGG", "ACAAAAAGAGAAGGG", "ATAGAAAAGGAAGGG", "AAGAAAAAAGAGAGG", "CAAGCTAAAGGGGGG"}
        // Consensus: AAAAAAAAAGAGGGG, Score: 39

        // --- Gibbs Sampler Search ---
        start = System.currentTimeMillis();
        motifs = gibbsSearch(DNA, 15, 200, 2000);
        end = System.currentTimeMillis();
        System.out.println(Arrays.toString(motifs));
        System.out.println("Score: " + Week3.score(motifs));
        System.out.println((end - start)/1000 + "s");
        System.out.println(Week3.bestMotif(motifs));

        // Best: [GAGAAAGGGAGGGCT, AAAAAATAGCAGGGT, AATAAACAGCGGGGT, CAGAAAAAAAGGGGT, ATAAAGTAGAGGGGG, AGAGGCGAGACGGGT, CAAAAAGAGAAGGGG, AAAAACCAGCGCGTG, AAAAAAGAGAGGAGT, AAGCTAAAGGGGGGT]
        // Consensus: AAAAAAGAGAGGGGT, Score: 38

    }

    public static String[] randomizedMotifSearch(String[] DNA, int k, int N) {
        String[] bestMotifs = randomSearch(DNA, k);

        for (int i = 0; i < N; i++) {
            String[] motifs = randomSearch(DNA, k);
            
            if (Week3.score(motifs) < Week3.score(bestMotifs)) {
                trace(Arrays.toString(motifs));
                trace(Arrays.toString(bestMotifs));
                trace("Better score found! " + Week3.score(motifs) + " < " + Week3.score(bestMotifs));
                bestMotifs = motifs.clone();
            }
        }
        return bestMotifs;
    }

    // Runs one 'instance' of motif search
    // (Run multiple times to get a good result)
    // Assumes all Strings in DNA are the same length
    private static String[] randomSearch(String[] DNA, int k) {
        int t = DNA.length;

        String[] motifs = new String[t];
        // Randomly select k-mers in each string from DNA
        for (int i = 0; i < t; i++) {
            int random = (int) (Math.random() * (DNA[i].length() - k));
            motifs[i] = DNA[i].substring(random, random+k);
        }
        String[] bestMotifs = motifs.clone();

        while (true) {
            // See Week 3 for explanation of profile and laplaceProfile
            double[][] profile = Week3.laplaceProfile(Week3.count(motifs), motifs.length);
            motifs = Week3.mostProbableArray(DNA, k, profile);

            if (Week3.score(motifs) < Week3.score(bestMotifs)) {
                bestMotifs = motifs.clone();
            } else {
                return motifs;
            }
        }

    }

    /*GibbsSampler(Dna, k, t, N)
    randomly select k-mers Motifs = (Motif1, …, Motift) in each string from Dna
    BestMotifs ← Motifs
    for j ← 1 to N
        i ← Random(t)
        Profile ← profile matrix constructed from all strings in Motifs except for Motifi
        Motifi ← Profile-randomly generated k-mer in the i-th sequence
        if Score(Motifs) < Score(BestMotifs)
            BestMotifs ← Motifs
    return BestMotifs */
    public static String[] gibbsSearch(String[] DNA, int k, int N, int repeat) {
        String[] bestMotifs = gibbsSampler(DNA, k, N);

        for (int i = 0; i < repeat; i++) {
            String[] motifs = gibbsSampler(DNA, k, N);
            
            if (Week3.score(motifs) < Week3.score(bestMotifs)) {
                // trace(Arrays.toString(motifs));
                // trace(Arrays.toString(bestMotifs));
                // trace("Better score found! " + Week3.score(motifs) + " < " + Week3.score(bestMotifs));
                bestMotifs = motifs.clone();
            }
        }
        return bestMotifs;
    }

    private static String[] gibbsSampler(String[] DNA, int k, int N) {
        int t = DNA.length;

        String[] motifs = new String[t];
        // Randomly select k-mers in each string from DNA
        for (int i = 0; i < t; i++) {
            int random = (int) (Math.random() * (DNA[i].length() - k));
            motifs[i] = DNA[i].substring(random, random+k);
        }
        
        String[] bestMotifs = motifs.clone();
        

        for (int count = 0; count < N; count++) {
            int i = (int) (Math.random() * t);
            String[] motifsNew = new String[t-1];
            for (int index = 0; index < t; index++) {
                if (index < i) {
                    motifsNew[index] = motifs[index];
                } else if (index > i) {
                    motifsNew[index-1] = motifs[index];
                }
            }
            double[][] profile = Week3.laplaceProfile(Week3.count(motifsNew), t);
            String generated = stringGenerator(DNA[i], k, profile);
            motifs[i] = generated;

            if (Week3.score(motifs) < Week3.score(bestMotifs)) {
                bestMotifs = motifs.clone();
            }
        }
        return bestMotifs;


    }

    public static String stringGenerator(String text, int k, double[][] profile) {
        double[] weights = new double[text.length() - k + 1];
        for (int position = 0; position <= text.length()-k; position++) {
            String kmer = text.substring(position, position+k);
            weights[position] = Week3.probabilty(kmer, profile);
        }
        int random = weightedRandom(weights);
        return text.substring(random, random+k);
        // String output = "";
        // for (int i = 0; i < profile[0].length; i++) {
        //     int random = weightedRandom(weights);
        //     output += BioFunctions.nucleotides[random];
        // }
        // return output;
    }

    public static int weightedRandom(double[] array) {
        double sum = 0;
        for (double d : array) {
            sum += d;
        }
        for (int i = 0; i < array.length; i++) {
            array[i] = array[i] / sum;
            
        }
        double random = Math.random();
        double counter = 0.0;
        for (int index = 0; index < array.length; index++) {
            counter += array[index];
            if (counter >= random) {
                return index;
            }
        }
        return 0;
    }

    private static void trace(Object object) {
        if (tracing) {
            System.out.println(object.toString());
        }
    }
}
