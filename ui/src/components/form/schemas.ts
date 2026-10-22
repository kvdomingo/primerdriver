import { z } from "zod";

export const settingsSchema = z.object({
  Tm_range_min: z.number(),
  Tm_range_max: z.number(),
  gc_range_min: z.number(),
  gc_range_max: z.number(),
  length_min: z.number().int().positive(),
  length_max: z.number().int().positive(),
  flank5_range_min: z.number().int().nonnegative(),
  flank5_range_max: z.number().int().nonnegative(),
  flank3_range_min: z.number().int().nonnegative(),
  flank3_range_max: z.number().int().nonnegative(),
  forward_overlap5: z.number().int().nonnegative(),
  forward_overlap3: z.number().int().nonnegative(),
  terminate_gc: z.boolean(),
  center_mutation: z.boolean(),
  primer_mode: z.enum(["complementary", "overlapping"]),
  expression_system: z.string().min(1),
});

export type Settings = z.infer<typeof settingsSchema>;

export const defaultSettings: Settings = {
  Tm_range_min: 75,
  Tm_range_max: 85,
  gc_range_min: 40,
  gc_range_max: 60,
  length_min: 25,
  length_max: 45,
  flank5_range_min: 11,
  flank5_range_max: 21,
  flank3_range_min: 11,
  flank3_range_max: 21,
  forward_overlap5: 9,
  forward_overlap3: 9,
  terminate_gc: true,
  center_mutation: true,
  primer_mode: "complementary",
  expression_system: "Homo sapiens",
};

const mutationType = z.enum(["sub", "ins", "del"]);

export const characterizeSchema = z.object({
  sequence: z.string().min(1),
  mismatched_bases: z.number().int().positive(),
  mutation_type: mutationType,
});

export const designSchema = z
  .object({
    sequence: z.string().min(1),
    mutation_type: mutationType,
    target: z.string(),
    position: z.number().int().nonnegative(),
    replacement: z.string(),
    settings: settingsSchema,
  })
  .refine(v => v.mutation_type === "ins" || v.target.length > 0, {
    path: ["target"],
    message: "Target is required",
  })
  .refine(v => v.mutation_type === "del" || v.replacement.length > 0, {
    path: ["replacement"],
    message: "Replacement is required",
  });
