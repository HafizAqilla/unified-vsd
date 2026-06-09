% extract_tables_text.m
try
    pdf_path = 'C:\Users\asus\Documents\SKRIPSI FILE\patient_specific_size_and_age_scaling_in_a_zero.776.pdf';
    txt_path = 'C:\Users\asus\Documents\VSD Main\unified-vsd-main\extracted_tables.txt';
    
    pdfFile = java.io.File(pdf_path);
    pdfDoc = org.apache.pdfbox.pdmodel.PDDocument.load(pdfFile);
    pdfStripper = org.apache.pdfbox.text.PDFTextStripper();
    
    full_txt = '';
    for page = 1:9
        try
            pdfStripper.setStartPage(page);
            pdfStripper.setEndPage(page);
            java_txt = pdfStripper.getText(pdfDoc);
            matlab_txt = char(java_txt);
            fprintf('Page %d extracted successfully (%d chars)\n', page, length(matlab_txt));
            full_txt = [full_txt, sprintf('\n\n--- PAGE %d ---\n', page), matlab_txt];
        catch page_err
            fprintf('Page %d extraction failed: %s\n', page, page_err.message);
        end
    end
    
    pdfDoc.close();
    
    % Fixed fopen syntax: 'n' specifies the native machine format, allowing 'UTF-8' encoding
    fid = fopen(txt_path, 'w', 'n', 'UTF-8');
    if fid == -1
        error('Failed to open output file for writing');
    end
    fprintf(fid, '%s', full_txt);
    fclose(fid);
    fprintf('Saved all successfully extracted pages to extracted_tables.txt\n');
    
catch ME
    fprintf('Global error: %s\n', ME.message);
    if exist('pdfDoc', 'var') && ~isempty(pdfDoc)
        pdfDoc.close();
    end
end
exit;
