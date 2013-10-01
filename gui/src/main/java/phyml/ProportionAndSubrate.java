package phyml;

import javax.swing.*;
import javax.swing.text.AttributeSet;
import javax.swing.text.BadLocationException;
import javax.swing.text.Document;
import javax.swing.text.PlainDocument;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

/**
 * Implements all components for specifying the "Proportion invariable sites"
 * and "One Substitutuon Rate Category" variable.
 *
 * @author Christoph Knapp
 *
 */

public class ProportionAndSubrate extends JPanel implements ActionListener {
	/**
	 * default id
	 */
	private static final long serialVersionUID = 1L;
	private CustomTextField parameter;
	private JRadioButton choice1;
	private JRadioButton choice2;
	private JRadioButton choice3;
	private JRadioButton choice4;

	/**
	 * Constructor method.
	 */
	public ProportionAndSubrate() {
		parameter = new CustomTextField("0.00");
		choice1 = new JRadioButton("Estimated");
		choice2 = new JRadioButton("Fixed");
		choice3 = new JRadioButton("No");
		choice4 = new JRadioButton("Yes");
		choice1.addActionListener(this);
		choice2.addActionListener(this);
		choice3.addActionListener(this);
		choice4.addActionListener(this);
		choice1.setSelected(true);

		CustomGridLayout layout = new CustomGridLayout();
		setLayout(layout);
		layout.setDimensions(1, 0.1);
		add(new JPanel());
		layout.setDimensions(0.01, 0.4);
		add(new JPanel());
		layout.setDimensions(0.33, 0.4);
		add(new JLabel("Proportion Invariable Sites"));
		layout.setDimensions(0.27, 0.4);
		add(parameter);
		layout.setDimensions(0.38, 0.4);
		JPanel p1 = new JPanel();
		CustomGridLayout lo1 = new CustomGridLayout();
		p1.setLayout(lo1);
		add(p1);
		lo1.setDimensions(0.2, 1);
		p1.add(new JPanel());
		lo1.setDimensions(0.4, 1);
		p1.add(choice1);
		lo1.setDimensions(0.4, 1);
		p1.add(choice2);
		layout.setDimensions(0.01, 0.4);
		add(new JPanel());
		layout.setDimensions(1, 0.1);
		add(new JPanel());
		layout.setDimensions(0.01, 0.4);
		add(new JPanel());
		layout.setDimensions(0.33, 0.4);
		add(new JLabel("Variable Rates Across Sites"));
		layout.setDimensions(0.27, 0.4);
		add(new JPanel());
		layout.setDimensions(0.38, 0.4);
		JPanel p2 = new JPanel();
		CustomGridLayout lo2 = new CustomGridLayout();
		p2.setLayout(lo2);
		add(p2);
		lo2.setDimensions(0.2, 1);
		p2.add(new JPanel());
		lo2.setDimensions(0.4, 1);
		p2.add(choice3);
		lo2.setDimensions(0.4, 1);
		p2.add(choice4);

        // setting variable rates to "Yes" by default
        choice4.setSelected(true);


	}

	/**
	 * CustomTextField extends JTextField. It is customized that the user is
	 * only able to type in Double values.
	 *
	 * @author Christoph Knapp
	 *
	 */
	class CustomTextField extends JTextField {
		/**
		 * default id
		 */
		private static final long serialVersionUID = 1L;

		/**
		 * Constructor Method.
		 *
		 * @param string
		 *            String : default value displayed in the text field.
		 */
		public CustomTextField(String string) {
			super(string);
			if (!isNumber(string)) {
				this.setText("");
			}
		}

		@Override
		protected Document createDefaultModel() {
			return new doubleOnlyDocument();
		}

		/**
		 * Tests wheter the input String is convertable to a Number.
		 *
		 * @param str
		 *
		 * @return
		 */
		private boolean isNumber(String str) {
			try {
				Double.parseDouble(str);
				return true;
			} catch (NumberFormatException e) {
				if (str.equals(".") && !this.getText().contains(".")) {
					return true;
				}
				return false;
			}
		}

		/**
		 * Customizes the default PainDocument to only accept numbers and one
		 * ".". Therefore only Double values can be typed in.
		 *
		 * @author Christoph Knapp
		 *
		 */
		class doubleOnlyDocument extends PlainDocument {
			/**
			 * default id
			 */
			private static final long serialVersionUID = 1L;

			@Override
			public void insertString(int offs, String str, AttributeSet a)
					throws BadLocationException {
				if (str == null || !isNumber(str)) {
					return;
				}
				super.insertString(offs, str, a);
			}
		}
	}

	@Override
	public void actionPerformed(ActionEvent e) {
		if (e.getSource() == choice1) {
			if(choice1.isSelected()){
				choice2.setSelected(false);
				parameter.setEnabled(false);
			}else{
				choice2.setSelected(true);
				parameter.setEnabled(true);
			}
		}
		if(e.getSource() == choice2){
			if(choice2.isSelected()){
				choice1.setSelected(false);
				parameter.setEnabled(true);
			}else{
				choice1.setSelected(true);
				parameter.setEnabled(false);
			}
		}
		if (e.getSource() == choice3) {
			if(choice3.isSelected()){
				choice4.setSelected(false);
				PhymlPanel.sM.setNumSubCatGammaAverageCompVisible(true);
			}else{
				choice4.setSelected(true);
				PhymlPanel.sM.setNumSubCatGammaAverageCompVisible(false);
			}
		}
		if(e.getSource() == choice4){
			if(choice4.isSelected()){
				choice3.setSelected(false);
				PhymlPanel.sM.setNumSubCatGammaAverageCompVisible(false);
			}else{
				PhymlPanel.sM.setNumSubCatGammaAverageCompVisible(true);
				choice3.setSelected(true);
			}
		}
	}

	/**
	 * Retrieves the value for the command line value for the proportion of
	 * invariable sites.
	 *
	 * @return String : "e" if estimated or Double value in String format if
	 *         fixed.
	 */
	public String getPropInvarSites() {
		if (choice1.isSelected()) {
			return "e";
		} else {
			return Double.parseDouble(parameter.getText()) + "";
		}
	}
}
