import sys
import re
from PyQt6.QtWidgets import (QApplication, QWidget, QLabel, QMessageBox,
                             QLineEdit, QTextEdit, QPushButton, QHBoxLayout, QVBoxLayout,
                             QFileDialog)
from PyQt6.QtCore import Qt


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
        self.gff_edit.setPlaceholderText("Select a GFF file...")
        self.gff_browse_button = QPushButton("Browse", self)
        self.gff_browse_button.clicked.connect(self.browse_gff_file)

        # Exon search fields
        search_label = QLabel("Search the position for an exon:")
        self.search_edit = QLineEdit()
        self.search_edit.setClearButtonEnabled(True)
        self.search_edit.setPlaceholderText("Provide an exon location")

        self.search_button = QPushButton("Search", self)
        self.search_button.clicked.connect(self.search_exon)
        self.search_edit.returnPressed.connect(self.search_exon)

        # Results fields
        results_label = QLabel("Information in TransDecoder:")
        results_label.setAlignment(Qt.AlignmentFlag.AlignLeft)
        self.chromosome_label = QLabel("Chr: ")
        self.strand_label = QLabel("Strand: ")
        self.results_text_edit = QTextEdit()
        self.results_text_edit.setReadOnly(True)
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
        gff_h_box.addWidget(self.gff_browse_button)

        search_h_box = QHBoxLayout()
        search_h_box.addWidget(self.search_edit)
        search_h_box.addWidget(self.search_button)

        info_h_box = QHBoxLayout()
        info_h_box.addWidget(self.chromosome_label)
        info_h_box.addWidget(self.strand_label)

        main_v_box = QVBoxLayout()
        main_v_box.addWidget(header_label)
        main_v_box.addLayout(gff_h_box)
        main_v_box.addWidget(search_label)
        main_v_box.addLayout(search_h_box)
        main_v_box.addWidget(results_label)
        main_v_box.addLayout(info_h_box)
        main_v_box.addWidget(self.results_text_edit)

        self.setLayout(main_v_box)

    # Slots for handling events
    def browse_gff_file(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Select GFF File", "", "GFF Files (*.gff *.gff3);;All Files (*)")
        if path:
            self.gff_edit.setText(path)

    def search_exon(self):
        file = self.gff_edit.text().strip()
        query_location = self.search_edit.text().strip()

        if not file or not query_location:
            QMessageBox.warning(self, "Input Required",
                                "<p>Please provide both a GFF file and a query location.</p>",
                                QMessageBox.StandardButton.Ok)
            return

        try:
            # Use mRNA ID as key so each isoform accumulates its own CDS list
            mrnas = {}

            with open(file) as f:
                for line in f:
                    if line.startswith("#") or line.strip() == "":
                        continue
                    cols = line.strip().split("\t")
                    if len(cols) < 9:
                        continue
                    feature = cols[2]
                    if feature == "mRNA":
                        id_match = re.search(r"ID=([^;]+)", cols[8])
                        if id_match:
                            mrna_id = id_match.group(1)
                            mrnas[mrna_id] = {
                                "chromosome": cols[0],
                                "strand": cols[6],
                                "isoform_id": mrna_id,
                                "cds_list": []
                            }
                    elif feature == "CDS":
                        parent_match = re.search(r"Parent=([^;]+)", cols[8])
                        if parent_match:
                            parent_id = parent_match.group(1)
                            if parent_id in mrnas:
                                mrnas[parent_id]["cds_list"].append((cols[3], cols[4]))

            matched = [iso for iso in mrnas.values()
                       if iso["cds_list"] and
                       any(query_location == s or query_location == e
                           for s, e in iso["cds_list"])]

            if not matched:
                self.chromosome_label.setText("Chr: ")
                self.strand_label.setText("Strand: ")
                self.results_text_edit.setText("Query Not Found...")
            else:
                output_parts = []
                for iso in matched:
                    label = iso["isoform_id"]
                    exon_str = ",".join(f"{s}..{e}" for s, e in iso["cds_list"])
                    output_parts.append(f"[{label}]\n{exon_str}")

                chromosomes = list(dict.fromkeys(iso["chromosome"] for iso in matched))
                strands = list(dict.fromkeys(iso["strand"] for iso in matched))
                self.chromosome_label.setText(f"Chr: {', '.join(chromosomes)}")
                self.strand_label.setText(f"Strand: {', '.join(strands)}")
                self.results_text_edit.setText("\n\n".join(output_parts))

        except FileNotFoundError as error:
            QMessageBox.warning(self, "Error",
                                f"<p>File not found.</p><p>{error}</p>",
                                QMessageBox.StandardButton.Ok)
        except OSError as error:
            QMessageBox.warning(self, "Error",
                                f"<p>Could not read file.</p><p>{error}</p>",
                                QMessageBox.StandardButton.Ok)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    sys.exit(app.exec())
