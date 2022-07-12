import java.io.IOException;

public class main {
    public static void main(String[] args) throws IOException {
        CHUI_PSO alg = new CHUI_PSO();
        alg.run();
        alg.printStats();
    }
}
