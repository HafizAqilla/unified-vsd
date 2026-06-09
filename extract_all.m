% extract_all.m
try
    pdf_path = 'C:\Users\asus\Documents\SKRIPSI FILE\patient_specific_size_and_age_scaling_in_a_zero.776.pdf';
    txt_path = 'C:\Users\asus\Documents\VSD Main\unified-vsd-main\pdf_extracted_text.txt';
    
    pdfFile = java.io.File(pdf_path);
    pdfDoc = org.apache.pdfbox.pdmodel.PDDocument.load(pdfFile);
    pdfStripper = org.apache.pdfbox.text.PDFTextStripper();
    
    java_txt = pdfStripper.getText(pdfDoc);
    matlab_txt = char(java_txt);
    
    pdfDoc.close();
    
    fid = fopen(txt_path, 'w', 'encoding', 'UTF-8');
    fprintf(fid, '%s', matlab_txt);
    fclose(fid);
    fprintf('Success! Full PDF text extracted (%d characters) to pdf_extracted_text.txt\n', length(matlab_txt));
catch ME
    fprintf('Error: %s\n', ME.message);
    if exist('pdfDoc', 'var') && ~isempty(pdfDoc)
        pdfDoc.close();
    end
end
exit;
