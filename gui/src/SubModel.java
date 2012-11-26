import java.awt.Color;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.FocusEvent;
import java.awt.event.FocusListener;
import java.io.File;

import javax.swing.BoxLayout;
import javax.swing.DefaultComboBoxModel;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.text.AttributeSet;
import javax.swing.text.BadLocationException;
import javax.swing.text.Document;
import javax.swing.text.PlainDocument;

/**
 * Class for instantiating all components to specify the actual substitution
 * model as well as setting there size and location.
 * 
 * @author Christoph Knapp
 * 
 */

public class SubModel extends JPanel implements ActionListener {
	/**
	 * default id
	 */
	private static final long serialVersionUID = 1L;
	private JComboBox<String> modelBox;
	private JLabel sMLab;
	private JLabel cMLab;
	private JLabel eFLab;
	private JLabel oRPLab;
	private CustomTextField curTextField;
	private JComboBox<String> equiBox;
	private JComboBox<String> optRateBox;
	private FreqPanel freqPanel;
	private String modelPath;
	private String molecularType;

	/**
	 * Constructor method for instantiating all components for specifying the
	 * substitution model and setting there size and location.
	 * 
	 * @param molecularType
	 *            Setting either "DNA" as molecule type or "AA".
	 */
	public SubModel(String molecularType) {
		this.molecularType = molecularType;
		modelPath = "";
		curTextField = new CustomTextField();
		String[] models;
		String[] equiChoices;
		if (molecularType.equals("DNA")) {
			models = new String[] { "HKY85", "F84", "TN93", "GTR", "custom",
					"JC69", "K80" };
			equiChoices = new String[] { "empirical", "optimised" ,"user defined" };
			eFLab = new JLabel("Nucleotide frequencies");
		} else {
			models = new String[] { "LG", "WAG", "Dayhoff", "JTT", "Blossum62",
					"Mt Rev", "Rt Rev", "Cp Rev", "DcMut", "VT", "Mt Mam",
					"Mt Art", "HIVw", "HIVb", "Read from file" };
			equiChoices = new String[] { "model", "empirical" };
			eFLab = new JLabel("Amino-acid frequencies");
		}
		modelBox = new JComboBox<String>(new DefaultComboBoxModel<String>(models));
		modelBox.addActionListener(this);
		equiBox = new JComboBox<String>(equiChoices);
		equiBox.addActionListener(this);
		String[] optRateChoices = new String[] { "yes", "no"};
		optRateBox = new JComboBox<String>(optRateChoices);
		setOptimiseRateOff(true);
		setOptimiseRate("NO");
		freqPanel = new FreqPanel();
    freqPanel.setCompVisible(false);
		sMLab = new JLabel("Substitution model");
		cMLab = new JLabel("Current model");
		oRPLab = new JLabel("Optimise rate parameter");
		if (molecularType.equals("AA")) {
			setCompVisible(false);
		}
		CustomGridLayout layout = new CustomGridLayout();
		setLayout(layout);
		layout.setDimensions(0.01, 1);
		add(new JPanel());
		layout.setDimensions(0.98, 0.04);
		add(new JPanel());
		layout.setDimensions(0.01, 1);
		add(new JPanel());
		layout.setDimensions(0.45, 0.29);
		JPanel p1 = new JPanel();
		add(p1);
		CustomGridLayout lo1 = new CustomGridLayout();
		p1.setLayout(lo1);
		lo1.setDimensions(0.7, 1);
		p1.add(sMLab);
		lo1.setDimensions(0.3, 1);
		p1.add(modelBox);
		layout.setDimensions(0.08, 0.29);
		add(new JPanel());
		layout.setDimensions(0.45, 0.29);
		JPanel p2 = new JPanel();
		add(p2);
		CustomGridLayout lo2 = new CustomGridLayout();
		p2.setLayout(lo2);
		lo2.setDimensions(0.72, 1);
		p2.add(cMLab);
		lo2.setDimensions(0.21, 1);
		p2.add(curTextField);
		layout.setDimensions(0.98, 0.04);
		add(new JPanel());
		layout.setDimensions(0.45, 0.29);
		JPanel p4 = new JPanel();
		add(p4);
		CustomGridLayout lo3 = new CustomGridLayout();
		p4.setLayout(lo3);
		lo3.setDimensions(0.7, 1);
		p4.add(eFLab);
		lo3.setDimensions(0.3, 1);
		p4.add(equiBox);
		layout.setDimensions(0.08, 0.29);
		add(new JPanel());
		layout.setDimensions(0.45, 0.29);
		JPanel p5 = new JPanel();
		add(p5);
		CustomGridLayout lo4 = new CustomGridLayout();
		p5.setLayout(lo4);
		lo4.setDimensions(1, 0.1);
		p5.add(new JPanel());
		lo4.setDimensions(0.72, 0.9);
		p5.add(oRPLab);
		lo4.setDimensions(0.21, 0.9);
		p5.add(optRateBox);
		layout.setDimensions(0.98, 0.04);
		add(new JPanel());
		layout.setDimensions(0.98, 0.3);
		add(freqPanel);
		setSelectedCompEnabled(2);
	}

	/**
	 * This method sets all components in the gui to visible or not visible.
	 * 
	 * @param b
	 *            boolean : true if visible, false if not visible.
	 */
	private void setCompVisible(boolean b) {
		cMLab.setVisible(b);
		oRPLab.setVisible(b);
		curTextField.setVisible(b);
		optRateBox.setVisible(b);
		//freqPanel.setCompVisible(b);
	}

	/**
	 * Changes the Comboboxmodel for the substitution model dropdown menu to DNA
	 * models or AA models.
	 * 
	 * @param molType
	 *            String : either "DNA" or "AA".
	 */
	public void setMoleculeType(String molType) {
		molecularType = molType;
		if (molecularType.equals("DNA")) {
			modelBox.setModel(new DefaultComboBoxModel<String>(new String[] { "HKY85",
					"F84", "TN93", "GTR", "custom", "JC69", "K80" }));
			setCompVisible(true);
			equiBox.setModel(new DefaultComboBoxModel<String>(new String[] { "empirical", 
					"optimised" ,"user defined" }));
			eFLab.setText("Nucleotide frequencies");
		} else if (molecularType.equals("AA")) {
			modelBox.setModel(new DefaultComboBoxModel<String>(new String[] { "LG",
					"WAG", "Dayhoff", "JTT", "Blossum62", "Mt Rev", "Rt Rev",
					"Cp Rev", "DcMut", "VT", "Mt Mam", "Mt Art", "HIVw",
					"HIVb", "Read from file" }));
			setCompVisible(false);
			equiBox.setModel(new DefaultComboBoxModel<String>(new String[] { "model",
					"empirical" }));
			eFLab.setText("Amino-acid frequencies");
		}
	}

	class CustomTextField extends JTextField implements ActionListener,
			FocusListener {
		/**
		 * default id
		 */
		private static final long serialVersionUID = 1L;

		public CustomTextField() {
			super("000000", 6);
			this.addActionListener(this);
			this.addFocusListener(this);
			this.setMinimumSize(new Dimension(40, 4));
			this.setMaximumSize(new Dimension(80, 25));
			this.setPreferredSize(new Dimension(60, 18));
		}

		protected Document createDefaultModel() {
			return new NumOnlyDocument();
		}

		/**
		 * 
		 * @author Christoph Knapp
		 * 
		 */
		class NumOnlyDocument extends PlainDocument {
			/**
			 * default id
			 */
			private static final long serialVersionUID = 1L;

			public void insertString(int offs, String str, AttributeSet a)
					throws BadLocationException {
				if (str == null || offs > 5 || !isNumber(str)
						|| this.getLength() > 5) {
					return;
				}
				char[] upper = str.toCharArray();
				super.insertString(offs, new String(upper), a);
			}

			private boolean isNumber(String str) {
				try {
					Integer.parseInt(str);
					return true;
				} catch (NumberFormatException e) {
					return false;
				}
			}
		}

		public void actionPerformed(ActionEvent e) {
			if (this.getText().length() < 6) {
				String txt = this.getText();
				for (int i = this.getText().length(); i < 6; i++) {
					txt = txt + "0";
				}
				this.setText(txt);
			}
		}

		public void focusGained(FocusEvent e) {
		}

		public void focusLost(FocusEvent e) {
			if (this.getText().length() < 6) {
				String txt = this.getText();
				for (int i = this.getText().length(); i < 6; i++) {
					txt = txt + "0";
				}
				this.setText(txt);
			}
		}
	}

	public void actionPerformed(ActionEvent e) {
		if (e.getSource() == modelBox) {
			if(modelBox.getSelectedItem().toString().equals("JC69")||
					modelBox.getSelectedItem().toString().equals("K80")){
				setOptimiseRateOff(true);
				setOptimiseRate("NO");
				PhymlPanel.setRatioBoxOff(false);
				equiBox.setSelectedIndex(2);
				setFrequencies("0.25","0.25","0.25","0.25");
				setSelectedCompEnabled(1);
				
			}else if (modelBox.getSelectedItem().toString().equals("custom")) {
				setOptimiseRateOff(true);
				setOptimiseRate("YES");
				PhymlPanel.setRatioBoxOff(false);
				setSelectedCompEnabled(true);
				setFrequencies("","","","");
			} else if (modelBox.getSelectedItem().toString().equals("Read from file")) {
				setOptimiseRateOff(true);
				setOptimiseRate("NO");
				JFileChooser fc = new JFileChooser();
				fc.setCurrentDirectory(new File(System.getProperty("user.dir")));
				int returnVal = fc.showOpenDialog(this);
				if (returnVal == JFileChooser.APPROVE_OPTION) {
					modelPath = fc.getSelectedFile().getAbsolutePath();
				}
				if (modelPath.equals("")) {
					modelBox.setSelectedIndex(0);
				}
				setFrequencies("","","","");
			}else if(modelBox.getSelectedItem().toString().equals("GTR")){
				setOptimiseRateOff(false);
				setOptimiseRate("NO");
				PhymlPanel.setRatioBoxOff(true);
				PhymlPanel.setRatioBoxDefault();
			}else {
				setOptimiseRateOff(true);
				setOptimiseRate("NO");
				PhymlPanel.setRatioBoxOff(false);
				modelPath = "";
				equiBox.setSelectedIndex(0);
				setSelectedCompEnabled(2);
				setFrequencies("","","","");
			}
		} else if (e.getSource() == equiBox) {
			if (equiBox.getSelectedItem().toString().equals("user defined")) {
				freqPanel.setCompEnabled(true);
				freqPanel.setCompVisible(true);
			} else {
				freqPanel.setCompEnabled(false);
				freqPanel.setCompVisible(false);
			}
		}
	}
	private void setOptimiseRate(String string) {
		if(string.equals("YES")){
			optRateBox.setSelectedIndex(0);
		}else if(string.equals("NO")){
			optRateBox.setSelectedIndex(1);
		}
	}

	private void setOptimiseRateOff(boolean b) {
		System.out.println("enabled = " + !b);
		optRateBox.setEnabled(!b);
	}

	/**
	 * Enables or disables components which should only be enabled if a user
	 * defined model is choosen.
	 * 
	 * @param i
	 *            int : which combination of enabled components to select.
	 */
	private void setSelectedCompEnabled(int i) {
		if(i==1){
			curTextField.setEnabled(false);
			equiBox.setEnabled(false);
			freqPanel.setCompEnabled(false);
			//freqPanel.setCompVisible(false);
		}else if(i==2){
			curTextField.setEnabled(false);
			equiBox.setEnabled(true);
			//freqPanel.setCompVisible(false);
      freqPanel.setCompEnabled(false);
		}
	}
	/**
	 * Enables or disables components which should only be enabled if a user
	 * defined model is choosen.
	 * 
	 * @param b
	 *            boolean : true if component enabled, false if disabled.
	 */
	private void setSelectedCompEnabled(boolean b) {
		curTextField.setEnabled(b);
		equiBox.setEnabled(b);
	}
	public void setFrequencies(String freqA, String freqC, String freqG, String freqT) {
		freqPanel.setFrequencies(freqA, freqC, freqG, freqT);
	}

	/**
	 * The selected model of substitution.
	 * 
	 * @return String : The selected mode i.e "HKY85". If custom is selected the
	 *         String represantation of the model consisting of 6 digits is
	 *         returned i.e. "000000".
	 */
	public String getSubModel() {
		if (molecularType.equals("DNA")
				&& modelBox.getSelectedItem().toString().equals("custom")) {
			return curTextField.getText();
		}
		return modelBox.getSelectedItem().toString();
	}

	/**
	 * If the moleeculetype AA is selected the user has the possibility to
	 * provide a rate file.
	 * 
	 * @return String : The full path to a specified Amino acid rate file.
	 */
	public String getAARateFileName() {
		return modelPath;
	}

	/**
	 * Retrieves the command line options for the specified model.
	 * 
	 * @return String : Comand line option of substitution model frequencies. -f
	 *         e, m, or “fA,fC,fG,fT” Nucleotide or amino-acid frequencies. – e
	 *         : the character frequencies are determined as follows : -
	 *         Nucleotide sequences: (Empirical) the equilibrium base
	 *         frequencies are estimated by counting the occurence of the
	 *         different bases in the alignment. - Amino-acid sequences:
	 *         (Empirical) the equilibrium amino-acid frequencies are estimated
	 *         by counting the occurence of the different amino-acids in the
	 *         alignment. – m : the character frequencies are determined as
	 *         follows : - Nucleotide sequences: (ML) the equilibrium base
	 *         frequencies are estimated using maximum likelihood. - Amino-acid
	 *         sequences: (Model) the equilibrium amino-acid frequencies are
	 *         estimated using the frequencies defined by the substitution
	 *         model. – “fA,fC,fG,fT” : only valid for nucleotide-based models.
	 *         fA, fC, fG and fT are floating numbers that correspond to the
	 *         frequencies of A, C, G and T respectively.
	 */
	public String getFrequencies() {
		if (molecularType.equals("DNA")) {
			if (equiBox.getSelectedItem().toString().equals("empirical")) {
				return "e";
			}else if(equiBox.getSelectedItem().toString().equals("optimised")){
				return "m";
			} else {
				return formatFrequencies(freqPanel.getFreq());
			}
		} else {
			if (equiBox.getSelectedItem().toString().equals("model")) {
				return "m";
			} else {
				return "e";
			}
		}
	}

	/**
	 * Formats the frequencies for the custom substitution model.
	 * 
	 * @param freq
	 *            String[] : array of frequencies.
	 * @return String : The frequencies in form of a comma seperated list for
	 *         passing them into the command line.
	 */
	private String formatFrequencies(String[] freq) {
		return freq[0] + "," + freq[1] + "," + freq[2] + "," + freq[3];
	}

	/**
	 * Whether the Rate parameter needs to be optimised.
	 * 
	 * @return boolean : true if the rate parameter should be optimised, false
	 *         otherwise.
	 */
	public boolean isOptimisedRateParameter() {
		if (optRateBox.getSelectedItem().toString().equals("yes")) {
			return true;
		}
		return false;
	}
}
