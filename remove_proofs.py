import re
import sys

def remove_proofs(tex_content):
    """Remove all proof environments from LaTeX content."""
    # Pattern to match \begin{proof}...\end{proof}
    pattern = r'\\begin\{proof\}.*?\\end\{proof\}'
    # Remove proof environments (including newlines)
    cleaned = re.sub(pattern, '', tex_content, flags=re.DOTALL)
    return cleaned

def main():
    if len(sys.argv) < 2:
        print("Usage: python remove_proofs.py <input.tex>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    
    # Read the tex file
    with open(input_file, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # Remove proofs
    cleaned_content = remove_proofs(content)
    
    # Generate output filename
    base_name = input_file.rsplit('.tex', 1)[0]
    output_file = f"{base_name}_no_proofs.tex"
    
    # Write the cleaned content
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(cleaned_content)
    
    print(f"Saved to {output_file}")

if __name__ == "__main__":
    main()