% inspect_pdf.m
try
    pdf_path = 'C:\Users\asus\Documents\SKRIPSI FILE\patient_specific_size_and_age_scaling_in_a_zero.776.pdf';
    pdfFile = java.io.File(pdf_path);
    pdfDoc = org.apache.pdfbox.pdmodel.PDDocument.load(pdfFile);
    
    fprintf('PDF loaded successfully. Page count: %d\n', pdfDoc.getNumberOfPages());
    
    fprintf('Trying to extract text...\n');
    pdfStripper = org.apache.pdfbox.text.PDFTextStripper();
    
    % Try to extract page 1 text only to see if that works
    pdfStripper.setStartPage(1);
    pdfStripper.setEndPage(1);
    txt = pdfStripper.getText(pdfDoc);
    fprintf('Page 1 extracted successfully! First 100 chars:\n%s\n', char(txt(1:min(100, length(txt)))));
    
    pdfDoc.close();
catch ME
    fprintf('Error Identifier: %s\n', ME.identifier);
    fprintf('Error Message: %s\n', ME.message);
    if isfield(ME, 'cause') && ~isempty(ME.cause)
        disp('Cause:');
        disp(ME.cause);
    end
    % If it is a Java exception, print Java stack trace
    if isa(ME, 'matlab.exception.JavaException')
        fprintf('Java Exception class: %s\n', class(ME.ExceptionObject));
        ME.ExceptionObject.printStackTrace();
    end
    for k = 1:numel(ME.stack)
        fprintf('Stack file: %s | line: %d | name: %s\n', ...
            ME.stack(k).file, ME.stack(k).line, ME.stack(k).name);
    end
    if exist('pdfDoc', 'var') && ~isempty(pdfDoc)
        pdfDoc.close();
    end
end
exit;
