import argparse
from cyvcf2 import VCF

def load_vcf_dict(filename):
    vcf_records = {}
    vcf = VCF(filename)
    for record in vcf:
        key = f"{record.CHROM}_{record.POS}_{record.REF}_{','.join(map(str, record.ALT))}"
        vcf_records[key] = record
    return vcf_records

def compare_vcf_records(record1, record2):
    differences = {}
    fields_to_compare = ['CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
    for field in fields_to_compare:
        value1 = getattr(record1, field, None)
        value2 = getattr(record2, field, None)
        if field == 'INFO':
            value1 = dict(record1.INFO)
            value2 = dict(record2.INFO)
            # Split PCGR_CSQ values for clearer comparison
            if 'PCGR_CSQ' in value1 or 'PCGR_CSQ' in value2:
                pcgr_csq1 = value1.get('PCGR_CSQ', '').split(',')
                pcgr_csq2 = value2.get('PCGR_CSQ', '').split(',')
                max_len = max(len(pcgr_csq1), len(pcgr_csq2))
                for i in range(max_len):
                    csq1 = pcgr_csq1[i] if i < len(pcgr_csq1) else ''
                    csq2 = pcgr_csq2[i] if i < len(pcgr_csq2) else ''
                    csq1_split = sorted(csq1.split('&'))
                    csq2_split = sorted(csq2.split('&'))
                    if csq1_split != csq2_split:
                        differences['PCGR_CSQ_order_diff'] = True
                        break
            # Compare other INFO fields individually
            for key in set(value1.keys()).union(set(value2.keys())):
                if key != 'PCGR_CSQ' and value1.get(key) != value2.get(key):
                    differences[key] = {'file1': value1.get(key, 'N/A'), 'file2': value2.get(key, 'N/A')}
        elif value1 != value2:
            differences[field] = {'file1': value1, 'file2': value2}
    return differences

def main():
    parser = argparse.ArgumentParser(description="Compare two VCF files record by record.")
    parser.add_argument("vcf_file1", help="Path to the first VCF file (original)")
    parser.add_argument("vcf_file2", help="Path to the second VCF file (new version)")
    args = parser.parse_args()

    vcf_file1 = args.vcf_file1
    vcf_file2 = args.vcf_file2

    records_vcf1 = load_vcf_dict(vcf_file1)
    records_vcf2 = load_vcf_dict(vcf_file2)

    records_only_in_vcf1 = set(records_vcf1.keys()) - set(records_vcf2.keys())
    records_only_in_vcf2 = set(records_vcf2.keys()) - set(records_vcf1.keys())
    common_records = set(records_vcf1.keys()).intersection(records_vcf2.keys())

    print(f"Records only in {vcf_file1} (original): {len(records_only_in_vcf1)}")
    print(f"Records only in {vcf_file2} (new version): {len(records_only_in_vcf2)}")
    print(f"Common records: {len(common_records)}")

    # Compare annotations for common records
    annotation_differences = {}
    order_diff_count_total = 0
    for key in common_records:
        record1 = records_vcf1[key]
        record2 = records_vcf2[key]
        differences = compare_vcf_records(record1, record2)
        if 'PCGR_CSQ_order_diff' in differences:
            order_diff_count_total += 1
            differences.pop('PCGR_CSQ_order_diff', None)
        if differences:
            annotation_differences[key] = differences

    # Output results
    print(f"Number of common records with different order annotations: {order_diff_count_total}")
    print(f"Number of common records with actual differences: {len(annotation_differences)}")

    # Output the count of records only in the new version
    print(f"Number of records only in the new version ({vcf_file2}): {len(records_only_in_vcf2)}")

if __name__ == "__main__":
    main()