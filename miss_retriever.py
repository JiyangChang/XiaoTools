import sys
import re
from PyQt6.QtWidgets import (QApplication, QWidget, QLabel, QMessageBox,
                             QLineEdit, QTextEdit, QPushButton, QHBoxLayout, QVBoxLayout)
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QFont


class MainWindow(QWidget):
    def __init__(self):
        super().__init__()
        self.initializeUI()

    def initializeUI(self):
        """Set up the application's GUI"""
        self.setGeometry(300, 300, 600, 300)
        self.setWindowTitle("XiaoTools - MR v0.1")

        self.setUpMainWindow()
        self.show()

    def setUpMainWindow(self):
        """Create and arrange widgets in the main window"""
        header_label = QLabel("Missing Retriever", self)
        # header_label.setFont(QFont("Arial", 16))
        header_label.setAlignment(Qt.AlignmentFlag.AlignCenter)

        # GFF fields
        gff_label = QLabel("GFF:")
        self.gff_edit = QLineEdit()
        self.gff_edit.setClearButtonEnabled(True)
        self.gff_edit.setPlaceholderText("Enter the path to the GFF file")

        # Exon search fields
        search_label = QLabel("Search the position for an exon:")
        self.search_edit = QLineEdit()
        self.search_edit.setClearButtonEnabled(True)
        self.search_edit.setPlaceholderText("Provide an exon location")

        self.search_button = QPushButton("Search", self)
        self.search_button.clicked.connect(self.search_exon)

        # Results fields
        results_label = QLabel("Information in TransDecoder:")
        results_label.setAlignment(Qt.AlignmentFlag.AlignLeft)
        self.chromosome_label = QLabel("Chr: ")
        self.strand_label = QLabel("Strand: ")
        self.results_text_edit = QTextEdit()
        self.results_text_edit.setPlaceholderText(
            """
            1. Run TransDecoder and get its result
            2. Copy the path to TransDecoder GFF file into GFF field
            3. Provide a boundary position of an missing exon to search field and click search
            4. Copy the result to GenomeView
            """
        )

        # Arrange layout of widgets
        gff_h_box = QHBoxLayout()
        gff_h_box.addWidget(gff_label)
        gff_h_box.addWidget(self.gff_edit)

        search_h_box = QHBoxLayout()
        search_h_box.addWidget(self.search_edit)
        search_h_box.addWidget(self.search_button)

        info_h_box = QHBoxLayout()
        info_h_box.addWidget(self.chromosome_label)
        info_h_box.addWidget(self.strand_label)

        results_h_box = QHBoxLayout()
        results_h_box.addWidget(self.results_text_edit)

        main_v_box = QVBoxLayout()
        main_v_box.addWidget(header_label)
        main_v_box.addLayout(gff_h_box)
        main_v_box.addWidget(search_label)
        main_v_box.addLayout(search_h_box)
        main_v_box.addWidget(results_label)
        main_v_box.addLayout(info_h_box)
        main_v_box.addLayout(results_h_box)

        self.setLayout(main_v_box)

    # Slots for handling events
    def search_exon(self):
        file = self.gff_edit.text()
        query_location = self.search_edit.text()
        cds_list = []

        try:
            with open(file) as f:
                for line in f:
                    if line != "\n":
                        line = line.strip().split("\t")
                        chromosome = line[0]
                        feature = line[2]
                        start_loc = line[3]
                        end_loc = line[4]
                        strand = line[6]
                        if feature == "CDS":
                            cds_list.append((start_loc, end_loc))
                    else:
                        final = []
                        found_flag = False
                        for i in cds_list:
                            if query_location in i:
                                found_flag = True
                        if found_flag:
                            for exon in cds_list:
                                start, end = exon
                                final.append(f"{start}..{end}")
                            self.strand_label.setText(f"Strand: {strand}")
                            self.chromosome_label.setText(f"Chr: {chromosome}")
                            self.results_text_edit.setText(",".join(final))
                            break
                        cds_list = []

        except FileNotFoundError as error:
            QMessageBox.warning(self, "Error",
                                f"""<p>File Not Found.</p><p>{error}</p>""",
                                QMessageBox.StandardButton.Ok)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    sys.exit(app.exec())
