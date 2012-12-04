import java.io.BufferedReader;
import java.io.InputStreamReader;

/**
 * Class CmdExec implements the call to phyml as a thread. This enables the
 * graphical user interface to get real time updates of the standard output.
 * 
 * @author Christoph Knapp
 * @version 27-June-2012
 * 
 */
public class CmdExec extends Thread {

	private String cmd;

	/**
	 * Constructor method for initialising an CmdExec object.
	 * 
	 * @param cmd
	 *            A command for calling phyml from terminal.
	 */
	public CmdExec(String cmd) {
		this.cmd = cmd;
	}

	@Override
	public void run() {
		int exitStatus = -1;
		try {
			Runtime rt = Runtime.getRuntime();
			Process process = rt.exec(cmd);
			BufferedReader input = new BufferedReader(new InputStreamReader(
					process.getInputStream()));
			String line1 = null;
			while ((line1 = input.readLine()) != null) {
				StandardOutPanel.setInput(line1);
                                Thread.sleep(2);
			}
			exitStatus = process.waitFor();
			System.out.println("Exit Status: " + exitStatus);
			PhymlPanel.SetSubmit(true);
			PhymlPanel.loadTrees();
		} catch (Throwable t) {
			t.printStackTrace();
		}
	}
}
