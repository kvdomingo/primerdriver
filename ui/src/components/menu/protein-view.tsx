import { useForm, useSelector } from "@tanstack/react-form";
import { $api } from "../../api";
import {
  AdvancedSettings,
  Form,
  MutationSelector,
  MutationType,
  ProteinSequenceInput,
  defaultSettings,
  designSchema,
} from "../form";
import Result from "./result";

function filterAminoAcids(value: string) {
  return value.toUpperCase().replace(/[^ARNDCQEGHILKMFPSTWYV]/g, "");
}

function ProteinView() {
  const mode = "PRO" as const;
  const { data: expressionSystems } = $api.useSuspenseQuery(
    "get",
    "/api/expressionsys",
  );
  const expressionList = expressionSystems.data;
  const primerDriver = $api.useMutation("post", "/api");

  const form = useForm({
    defaultValues: {
      sequence: "",
      mutation_type: "",
      target: "",
      position: 0,
      replacement: "",
      settings: defaultSettings,
    },
    validators: { onSubmit: designSchema },
    onSubmit: ({ value }) => primerDriver.mutate({ body: { mode, ...value } }),
  });
  const values = useSelector(form.store, state => state.values);
  const isValid = useSelector(
    form.store,
    state => designSchema.safeParse(state.values).success,
  );

  return (
    <>
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
      <Result mode={mode} status={primerDriver.status} data={primerDriver.data} />
    </>
  );
}

export default ProteinView;
