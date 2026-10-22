import { useForm, useSelector } from "@tanstack/react-form";
import { useNavigate } from "@tanstack/react-router";
import { useEffect, useState } from "react";
import api from "../../api";
import { type ResultData, usePrimerDriverStore } from "../../store";
import {
  AdvancedSettings,
  Form,
  MutationSelector,
  MutationType,
  ProteinSequenceInput,
  defaultSettings,
} from "../form";
import LoadingScreen from "../shared/LoadingScreen";

function filterAminoAcids(value: string) {
  return value.toUpperCase().replace(/[^ARNDCQEGHILKMFPSTWYV]/g, "");
}

function ProteinView() {
  const [expressionList, setExpressionList] = useState<string[]>([]);
  const navigate = useNavigate();
  const { setMode, setLoadedResultsFromModule, setResults } = usePrimerDriverStore();
  const mode = "PRO" as const;

  const form = useForm({
    defaultValues: {
      sequence: "",
      mutation_type: "",
      target: "",
      position: 0,
      replacement: "",
      settings: defaultSettings,
    },
    onSubmit: async ({ value }) => {
      let data: ResultData;
      try {
        data = (await api.data.primerDriver({ mode, ...value })).data;
      } catch (err) {
        console.error(err);
        data = "Request failed. Please try again later.";
      }
      setLoadedResultsFromModule(true);
      setResults({ data, loaded: true });
      await navigate({ to: "/results" });
    },
  });
  const isSubmitting = useSelector(form.store, state => state.isSubmitting);
  const values = useSelector(form.store, state => state.values);
  const isValid =
    values.sequence.length > 0 &&
    values.mutation_type !== "" &&
    (values.mutation_type === "sub"
      ? values.target.length > 0 && values.replacement.length > 0
      : values.target.length > 0 || values.replacement.length > 0);

  useEffect(() => {
    api.data
      .expressionSystem()
      .then(res => setExpressionList(res.data.data))
      .catch(err => console.error(err));

    setMode(mode);
  }, [setMode, mode]);

  if (isSubmitting) return <LoadingScreen />;
  return (
    <Form
      title="Protein-based Primer Design"
      onSubmit={form.handleSubmit}
      isValid={isValid}
    >
      <form.Field name="sequence">
        {field => (
          <ProteinSequenceInput
            sequence={field.state.value}
            onChange={value => field.handleChange(filterAminoAcids(value))}
          />
        )}
      </form.Field>
      <form.Field name="mutation_type">
        {field => (
          <MutationType
            handleChange={field.handleChange}
            mutation_type={field.state.value}
          />
        )}
      </form.Field>
      <MutationSelector
        handleChangeTarget={value =>
          form.setFieldValue("target", filterAminoAcids(value))
        }
        handleChangeReplacement={value =>
          form.setFieldValue("replacement", filterAminoAcids(value))
        }
        handleChangePosition={value => form.setFieldValue("position", value)}
        mutation_type={values.mutation_type}
        position={values.position}
        replacement={values.replacement}
        target={values.target}
      />
      <form.Field name="settings">
        {field => (
          <AdvancedSettings
            value={field.state.value}
            onChange={field.handleChange}
            mode={mode}
            expressionList={expressionList}
          />
        )}
      </form.Field>
    </Form>
  );
}

export default ProteinView;
