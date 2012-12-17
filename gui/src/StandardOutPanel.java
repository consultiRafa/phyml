import java.awt.Font;
import java.awt.Rectangle;

import javax.swing.DropMode;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextPane;

/**
 * Extends JPanel and implements the components to display the standard output
 * of phyml.
 * 
 * @author Christoph Knapp
 * 
 */

public class StandardOutPanel extends JPanel {
	/**
	 * default id
	 */
	private static final long serialVersionUID = 1L;
	private static JTextArea editorPane;

	/**
	 * Constructor method implements the components to display the standard
	 * output of phyml.
	 */
	public StandardOutPanel() {
		CustomGridLayout layout = new CustomGridLayout();
		setLayout(layout);
		layout.setDimensions(1, 1);
		editorPane = new JTextArea();
    DefaultCaret caret = (DefaultCaret)editorPane.getCaret();
		caret.setUpdatePolicy(DefaultCaret.ALWAYS_UPDATE);
		Font font = new Font("Courier", Font.PLAIN, 12);
		editorPane.setFont(font);
		editorPane.setEditable(false);
		editorPane.setDropMode(DropMode.INSERT);
		JScrollPane editorScrollPane = new JScrollPane(editorPane);
		editorScrollPane
				.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);
		editorScrollPane
				.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
		add(editorScrollPane);
	}

	/**
	 * Adds a new line to the standard output panel.
	 * 
	 * @param line
	 *            String : A line to add to the output panel
	 */
	public static void setInput(String line) {
		if(line.contains("\t")){
			System.out.println(line);
		}
		editorPane.append(line + "\n");
	}
}
