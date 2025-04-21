# Python script to predict eye color based on 8-plex SNP system
def predict_eye_color(genotypes):
    rs12913832 = genotypes.get('rs12913832', '')
    rs16891982 = genotypes.get('rs16891982', '')
    rs12203592 = genotypes.get('rs12203592', '')
    rs12896399 = genotypes.get('rs12896399', '')
    rs6119471  = genotypes.get('rs6119471', '')

    # Step 1
    if rs12913832 in ['AA', 'AG', 'GA']:
        if rs6119471 == 'GG' or rs16891982 == 'CC':
            return "Brown"
        elif rs12913832 == 'GA' and rs12203592 == 'TT':
            return "Green"
        else:
            step1 = "Not blue (Brown or Green)"
    elif rs12913832 == 'GG':
        if rs16891982 == 'CC':
            return "Green"
        elif rs12203592 == 'TT':
            return "Blue"
        else:
            step1 = "Not brown (Blue or Green)"
    else:
        step1 = "Inconclusive"

    # Step 2
    if step1 == "Not brown (Blue or Green)" and rs12913832 == 'GG' and rs12896399 == 'TT':
        return "Blue"
    elif step1 == "Not blue (Brown or Green)" and rs12913832 == 'AA' and rs12896399 == 'GG':
        return "Brown"

    return step1

# Example usage:
your_genotypes = {
    'rs12913832': 'AG',
    'rs16891982': 'CG',
    'rs12203592': 'CT',
    'rs12896399': 'GG',
    # 'rs6119471': '??' # Not available
}

predicted_eye_color = predict_eye_color(your_genotypes)
print("Predicted eye color:", predicted_eye_color)
