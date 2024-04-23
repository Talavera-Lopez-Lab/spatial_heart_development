import cellxgene_census

def get_value_filters(keys=None, organism='homo_sapiens'):
    with cellxgene_census.open_soma(census_version="2023-12-15") as census:
        if keys == None:
            keys = list(census['census_data'][organism].obs.keys())
        values_list = []
        for key in keys:
            values_df = (
                census['census_data'][organism].obs
                .read(column_names=[key])
                .concat()
                .to_pandas()
                .drop_duplicates()
                .to_dict('list')
            )
            values_list.append(values_df)
    values_dict = {k: v for dict in values_list for k, v in dict.items()}
    return values_dict